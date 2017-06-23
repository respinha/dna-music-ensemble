
from __future__ import division

import subprocess

import time
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio import motifs
from Bio.Align import AlignInfo, MultipleSeqAlignment

from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import MuscleCommandline


from scipy.stats import itemfreq

import numpy as np
import pandas

import os
import sys
import random

from Bio.Phylo.TreeConstruction import DistanceCalculator

from matplotlib import pyplot as plt

from music21 import note, scale, instrument, stream, duration, tempo
from music21.instrument import StringInstrument, WoodwindInstrument, BrassInstrument, PitchedPercussion, Instrument
from music21 import Music21Object

from algorithms import *
from config import GLOBALS, MIN_TEMPO


def read_phylo_tree(sequence_file):
    tree = Phylo.read(sequence_file.split('.')[0] + '.dnd', 'newkick')
    Phylo.draw_ascii(tree)


# returns array of possible motifs
def gen_motifs(fasta_file, **args):
    cmd = []
    cmd.append('meme')
    cmd.append(fasta_file)
    cmd.append('-dna')
    cmd.append('-mod')
    cmd.append('zoops')

    valid_args = ['nmotifs',  # max number of motifs
                  'nsites',  # number of sites for each motif
                  'prior',  # prior distribution file
                  'minw',  # minimum motif width
                  'maxw']  # maximum motif width

    output = 'meme_out'
    for key, value in args.iteritems():
        if key not in valid_args:
            raise Exception('Invalid argument: ', key)

        if key == 'o' or key == 'oc':
            output = value

        cmd.append("-" + key)
        cmd.append(str(value))

    with open(fasta_file, 'rU') as handle:

        subprocess.call(cmd)

        with open(output + '/meme.txt', 'rU') as meme_file:
            record = motifs.parse(meme_file, 'meme')
            return record


# aux function
# creates a generator from an iterable (for example, another generator)
# from the original iterable's first n elements
def generator_from_iterable(iterable, n):

    i = 0
    for r in iterable:
        if i < n:
            yield r
        else: break
        i += 1


# generates a MSA from a file with a set of sequences
# arguments can be:
#   seq_vector: vector specifying subset of sequences by reference
#   n_sequences: first n sequences of file
#   MSA algorithm (default: Clustal)
def gen_alignment(input_file, seq_vector=None, n_sequences=None, algorithm='mafft', output_file='output'):

    assert input_file is not None and os.path.isfile(input_file)
    assert output_file is not None

    assert seq_vector is not None or n_sequences is not None, \
        'Both arguments are None (sequence vector and number of sequences)'

    assert isinstance(seq_vector, list) or isinstance(n_sequences, int), \
        'Either one of two must be provided: sequence vector or number of sequences'

    assert algorithm in GLOBALS['SUPPORTED ALGORITHMS'], \
        'Algorithm does not match any of the currently supported MSA algorithms'

    assert isinstance(input_file, str)

    iterable = SeqIO.parse(open(input_file, 'rU'), 'fasta')

    tmp_file = 'pre_alignment.fna'

    if seq_vector is not None:
        sequences = (r for r in iterable if r.description.split('|')[-1] in seq_vector)
    else:
        sequences = generator_from_iterable(iterable, n_sequences)

    sequences = [x for x in sequences]
    if len(sequences) == 0:
        print 'No sequences were found'
        sys.exit(0)

    #print sequences
    SeqIO.write(sequences, tmp_file, 'fasta')

    try:

        t0 = time.time()
        if algorithm == 'clustal':

            if not output_file.endswith('.aln'):
                output_file += '.aln'

            algorithm = 'clustalw2'
            cline = ClustalwCommandline(algorithm,
                                        infile=tmp_file,
                                        outfile=output_file + '.aln')
        elif algorithm == 'muscle':

            if not output_file.endswith('.fna'):
                output_file += '.fna'

            alg = r"/usr/local/bin/muscle3.8.31_i86linux64"
            cline = MuscleCommandline(alg, input=tmp_file,
                                      out='source_sequences/muscle_' + str(n_sequences) + '.fna',
                                      clwstrict=True)
        elif algorithm == 'mafft':

            if not output_file.endswith('.fasta'):
                output_file += '.fasta'

            alg = r"/usr/local/bin/mafft"
            cline = MafftCommandline(alg,
                                     input=tmp_file,
                                     clustalout=True)
        else:
            print 'Unknown algorithm\n'
            sys.exit(0)

        stdout, stderr = cline()

        if algorithm == 'mafft':
            with open(output_file, "wb") as handle:
                handle.write(stdout)

        print 'Elapsed time: ' + str((time.time() - t0) / 60)

        return output_file
    except:
        print 'Error aligning with ' + algorithm
        

# retrieves a distance matrix from:
#   a) a multiple sequence alignment
#   b) a file containing a multiple sequence alignment
def get_distance_matrix(msa):

    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(msa)

    return np.array([row for row in distance_matrix])


# clusters all sequences in a MSA
def get_clusters_from_alignment(msa, **kwargs):
    
    assert isinstance(msa, MultipleSeqAlignment)

    print 'Retrieving distance matrix'
    dm = get_distance_matrix(msa)

    instruments = np.array(len(msa))

    if dm is not None:
        
        assert 'algorithm' in kwargs.keys(), 'No algorithm specified for clustering'
        
        algorithm = kwargs['algorithm']

        if 'nclusters' not in kwargs.keys():
            nclusters = len(msa) / 2
        else:
            nclusters = kwargs['nclusters']
            assert isinstance(nclusters, int)

        if algorithm == 'kmeans':
                
            from sklearn.cluster import KMeans
            
            model = KMeans(n_clusters=nclusters, random_state=0)
            model.fit(dm)

            clusters = model.labels_
            # centroids = model.cluster_centers_

        elif algorithm == 'hierarchical':

            from scipy.cluster.hierarchy import dendrogram, linkage
            from scipy.cluster.hierarchy import fcluster

            print 'Retrieving cluster tree'
    
            Z = linkage(dm)

            """if 'dendrogram' in kwargs.keys():
                if kwargs['dendrogram']:
    
                    plt.title('Hierarchical Clustering Dendrogram')
                    plt.xlabel('sample index')
                    plt.ylabel('distance')
        
                    dendrogram(
                        Z,
                        leaf_rotation=90.,  # rotates the x axis labels
                        leaf_font_size=8.,  # font size for the x axis labels
                    )
                    plt.savefig('dendogram.png')
            """

            max_d = 0.01 # TODO: make this dynamic
            clusters = fcluster(Z, max_d, criterion='distance')

        else:
            print 'Invalid cluster algorithm'
            raise NotImplementedError

        from random import shuffle

        instruments_pool = PitchedPercussion.__subclasses__()
        shuffle(instruments_pool)
        instruments_pool = instruments_pool[:nclusters]

        # TODO: TEST!
        i = 0
        for cluster in clusters:
            instruments[i] = instruments_pool[cluster]
            i += 1

        return instruments

    return None


# clusters all sequences in a MSA file
def cluster_alignment(alignment, depth=1):

    assert isinstance(alignment, MultipleSeqAlignment)

    print 'Retrieving clusters...'

    clusters = get_clusters_from_alignment(alignment, depth)

    n_clusters = max(clusters)
    if n_clusters <= 4:

        # numpy array containing pointers to each instrument family
        # and the respectively assigned instrument
        sequence_instruments = np.zeros(len(clusters), dtype=('uint8,uint8'))

        import random

        for i in range(0, len(clusters)):

            idx = clusters[i]
            family = GLOBALS['FAMILIES'][idx]
            print family

            try:
                instruments = family.__subclasses__()
            except TypeError:
                instruments = family.__subclasses__(family)

            rnd = random.randint(0, len(instruments) - 1)
            print instruments[rnd]  # TODO: build list of instruments by priorities

            sequence_instruments[i] = (idx, rnd)

    return get_clusters_from_alignment(alignment, depth)


# aux function
# converts a MSA file in any format to 'phylip-relaxed'
def msa_to_phylip(msa):
    assert os.path.isfile(msa), "MSA file does not exist: " + msa

    out_file = msa.split('.')[0] + '.phy'
    AlignIO.convert(msa, 'clustal', out_file, 'phylip-relaxed')

    return out_file


def gen_dynamics_vector(msa, dynamics_algorithm):

    # criteria: local, avg ou median entropy
    assert isinstance(dynamics_algorithm, DynamicsAlgorithm)
    assert 'window_size' in dynamics_algorithm.keys(), 'Empty window for dynamics algorithm'
    assert dynamics_algorithm['algorithm'] == DynamicsAlgorithm.SHANNON_INDEX # todo: only one option for now; implement simpson index afterwards

    window = dynamics_algorithm['window_size']

    if 'gap_threshold' not in dynamics_algorithm.keys():
        gap_threshold = 0.7
    else:
        gap_threshold = dynamics_algorithm['gap_threshold']

    if 'criteria' not in dynamics_algorithm.keys():
        criteria = 'local'
    else:
        criteria = dynamics_algorithm['criteria']

    if 'levels' not in dynamics_algorithm.keys():
        levels = 5
    else:
        levels = dynamics_algorithm['levels']

    aln_len = np.alen(msa[0])

    from math import ceil

    n_windows = ceil(float(aln_len)/window)+1
    gaps_below_thresh = np.zeros(n_windows, dtype=np.bool)

    window_idx = 0

    # first iterating through MSA to identify
    # which windows have a percentage of gaps
    # below the established threshold

    for i in range(0, aln_len, window):

        if i + window > aln_len: window = aln_len - i

        local_is_below_threshold = np.zeros((window,), dtype=np.bool)
        for j in range(i, i+window):

            column = np.array([c for c in msa[:,j]])

            n_gaps = np.count_nonzero(column == '-')

            if float(n_gaps) / len(column) < gap_threshold:
                local_is_below_threshold[j-i] = True

        n_ungapped_regions = np.alen(np.where(local_is_below_threshold)[0])

        if n_ungapped_regions < gap_threshold * window:
            gaps_below_thresh[window_idx] = True

        window_idx += 1

    from scipy.stats import entropy

    # dynamics_vector = np.zeros((n_windows,),
    dynamics_vector = np.zeros((n_windows,),
                               dtype=[('entropy', np.float), ('vol', np.float)])

    entropies_idx = 0

    for i in range(0, aln_len, window):

        if i + window > aln_len: window = aln_len - i

        # if this window has a percentage of gaps
        # above the considered threshold
        if not gaps_below_thresh[entropies_idx]:
            dynamics_vector['entropy'][entropies_idx] = -1
            continue

        local_entropy = np.zeros((window,))

        for j in range(i, i+window):

            column = msa[:, j]

            column_symbols = column[np.where(column != '-')]

            if len(column_symbols) < (1-gap_threshold) * len(column):

                # apply thresholding +
                pass
            else:
                # Shannon entropy of all column symbols

                counts = itemfreq(column_symbols)
                counts = np.array([float(count) for count in counts[:,1]])

                counts /= np.sum(counts)

                local_entropy[j-i] = entropy(counts, base=2)

        if criteria == 'local':
            dynamics_vector['entropy'][entropies_idx] = np.sum(local_entropy)
        elif criteria == 'average':
            dynamics_vector['entropy'][entropies_idx] = np.average(local_entropy)
        elif criteria == 'median':
            dynamics_vector['entropy'][entropies_idx] = np.median(local_entropy)
        else:
            print 'Unsupported criteria ' + str(criteria) + ' for entropy aggregation'
            sys.exit(1)

        entropies_idx += 1

    max_vol = 0.95
    min_vol = 0.30

    entropies = dynamics_vector['entropy']

    split_info = np.array_split(np.sort(np.unique(entropies)), levels)   # splitting info into classes
    volumes = np.linspace(min_vol, max_vol, num=levels)                # vector with all possible volumes

    for i in range(0, int(n_windows)):

        for j in range(0, np.alen(split_info)):
            if entropies[i] == -1:
                dynamics_vector['vol'][i] = -1
            elif entropies[i] <= split_info[j][-1]:

                dynamics_vector['vol'][i] = volumes[j]
                break

    return dynamics_vector


def add_dynamics_to_score(dynamics_vector, score, window_size, instruments, max_rest_tempo=3):

    # consistency check
    assert isinstance(score, stream.Score) and isinstance(dynamics_vector, np.ndarray)
    assert len(instruments) == len(score.parts)
    assert all(isinstance(i, instrument.Instrument) for i in instruments)

    vol_idx = 0
    score_tempo = score.getElementsByClass(tempo.MetronomeMark)[0]

    # creating filtered score
    # starting with empty parts
    final_score = stream.Score()
    length = np.inf

    for p in range(0, len(score.parts)):

        if len(score.parts[p]) < length:
            length = len(score.parts[p])

        part = stream.Part()
        part.insert(0, score_tempo)
        part.insert(0, instruments[p])

        final_score.append(part)

    # iterating through a music score in chunks of 'window_size' length

    for i in range(0, length, window_size):
        window_size = length - i if i + window_size > length else window_size

        print 'len', length, 'window_size', window_size

        if dynamics_vector[vol_idx] <= 0:
            r = note.Rest()

            for part in final_score.parts:
                part.append(r)
                part[-1].seconds = max_rest_tempo

        else:
            # iterating over parts
            for j in range(0, len(final_score.parts)):

                # scores have same number of parts
                final_part = final_score.parts[j]
                part = score.parts[j]

                for k in range(i, i + window_size):

                    n = part[k]
                    if isinstance(n, note.GeneralNote):  # if Note or Rest

                        final_part.append(n)
                        if isinstance(n, note.Note):  # if Note

                            final_part[-1].volume = dynamics_vector[vol_idx]

        vol_idx += 1

    print 'SECS', final_score[0].seconds

    assert len(final_score.parts) == len(score.parts)
    return final_score


# returns a tuple containing:
#   - a word distance vector per word (A,C,G,T)
#   - a label array with the assigned duration of each nucleotide
def gen_pitch_duration_vectors(sequence, pitch_algorithm, durations_algorithm):

    # kwargs:
    #   - window_duration
    #   - step
    #   - duration_mapping

    assert isinstance(pitch_algorithm, PitchAlgorithm) and isinstance(durations_algorithm, DurationsAlgorithm)
    assert sequence is not None

    # step is the number of columns/characters that are mapped
    length = np.alen(sequence)

    # auxiliary structures
    distance_vectors = dict()  # keys are nucleotides; values are np arrays
    last_occurrence = dict()  # aux dict used for word distances

    step = 1
    window = 1500
    window_duration = 8

    for key in durations_algorithm.keys():
        if key == 'n_nucleotides':
            step = durations_algorithm['n_nucleotides']
            assert isinstance(step, int) and step > 0

        elif key == 'window_size':
            window = durations_algorithm['window_size']
            assert isinstance(window, int) and window > 0

        elif key == 'window_duration':
            window_duration = durations_algorithm['window_duration']
            assert (isinstance(window_duration, float) or isinstance(window_duration, int)) and window > 0

    if 'n_nucleotides' in pitch_algorithm.keys():
        # TODO: implementar para N nucleotidos
        assert step == pitch_algorithm['n_nucleotides']

    d_algorithm = durations_algorithm['algorithm']
    p_algorithm = pitch_algorithm['algorithm']

    # durations vector
    durations = np.zeros((length,), dtype=np.float)

    for i in range(0, length, window):

        if i + window <= length:
            boundary = i + window
        else:
            boundary = i + (length - i)

        # retrieving region to calculate relative frequencies
        subset = sequence[i:boundary]

        # grouping frequencies of A,C,G,T  within a block
        counts = dict(itemfreq(subset))

        # This algorithm assigns durations dynamically to nucleotides
        # based on their relative frequency
        if d_algorithm == durations_algorithm.FREQUENCIES_DYNAMIC:

            counts_sum = 0

            for j in range(i, boundary):
                counts_sum += int(counts[subset[j-i]])

        total_time = 0.0

        for j in range(i, boundary, step):

            # pitch algorithm
            if p_algorithm == PitchAlgorithm.WORD_DISTANCES:

                # TODO: fazer uma primeira passagem pela sequencia inteira
                # para saber o numero exato de nucleotidos e reservar np.arrays??
                letter = sequence[j]

                if letter not in distance_vectors.keys():
                    distance_vectors[letter] = []

                else:
                    diff = j - last_occurrence[letter]
                    distance_vectors[letter].append(diff)

                """if d_mapping == 3:
                    idx = durations[letter][1]
                    durations[letter][0][idx] = diff

                    idx += 1
                    durations[letter][1] = idx

                    total_time += float(diff)"""

                last_occurrence[letter] = j

            # durations algorithm

            local_count = int(counts[subset[j - i]])
            freq = float(local_count) / len(subset)

            """
                TODO: considerar isto?
                if d_mapping == 0:

                duration_labels = {0.25: 'eighth', 0.5: 'quarter', 0.75: 'half', 1.0: 'whole'}
                local_duration = 'eighth'

                keys = duration_labels.keys()

                for k in keys:

                    if k > freq:
                        local_duration = duration_labels[k]
                        break
            """

            # frequency-biased algorithms
            if d_algorithm == DurationsAlgorithm.FREQUENCIES_DYNAMIC:

                assert counts_sum is not None and counts_sum > 0, 'Inconsistent values for counts of characters in region'

                local_duration = float(window_duration) * float(local_count) / float(counts_sum)
                total_time += local_duration

                # assert local_duration > MIN_TEMPO, \
                #    'Higher tempo required for each subsequence; too short duration was calculated: ' + str(local_duration)

                # todo: test
                if local_duration < MIN_TEMPO:
                    # print 'WARNING: Converting local duration to minimum tempo!'
                    local_duration = MIN_TEMPO

            # duration biased with discrete attributions
            elif d_algorithm == DurationsAlgorithm.FREQUENCIES_DISCRETE:

                duration_labels = np.array(['32nd','16th','eighth','quarter','half','whole'])
                frequencies = np.linspace(0.0, 0.5, num=len(duration_labels)-1)

                # dictionary mapping a discrete value for word frequencies to each duration label
                duration_labels = {frequencies[i]:duration_labels[i] for i in range(0, len(duration_labels)-1)}
                duration_labels['whole'] = 1.0

                local_duration = '32nd'
                keys = duration_labels.keys()

                # print freq
                for k in keys:
                    if k > freq:
                        local_duration = duration_labels[k]
                        break

                local_duration = duration.Duration(local_duration).quarterLength
                total_time += local_duration

            else:
                print 'Invalid mapping introduced: ', d_algorithm['algorithm']
                raise NotImplementedError

            durations[j] = local_duration

        # if d_mapping > 0 and d_mapping < 3:
        #     total = np.sum([durations[j] for j in range(i, boundary)])

        if d_algorithm == DurationsAlgorithm.FREQUENCIES_DISCRETE:

            if round(float(total_time), 5) != round(float(window_duration), 5):

                ratio = float(window_duration) / float(total_time)
                total_time = 0.0
                for j in range(i, boundary):
                    durations[j] *= ratio

                    assert durations[j] >= MIN_TEMPO,\
                        'Higher tempo required for each subsequence; too short duration was calculated: ' + str(durations[j])

                    total_time += durations[j]

        """elif d_mapping == 3:

            if round(float(total_time), 5) != round(float(window_duration), 5):

                ratio = float(window_duration) / float(total_time)
                total_time = 0.0

                print subset
                for key, array in iter(durations.items()):

                    print key
                    sub_array = array[0][i:boundary]

                    for j in range(i, boundary):

                        print sub_array[j]
                        sub_array[j] *= ratio
                        total_time += sub_array[j]

                        assert sub_array[j] >= MIN_TEMPO, \
                        'Higher tempo required for each subsequence; too short duration was calculated: ' + str(sub_array[j])


        if d_mapping > 0:
            assert round(float(total_time), 5) == round(float(window_duration), 5), \
                'Durations don\'t  match: ' + str(float(total_time)) + ' ' + str(float(window_duration))"""

    for x, y in iter(distance_vectors.items()):

        diff = length - last_occurrence[x]
        distance_vectors[x].append(diff)

    d_vectors_len = 0
    for x in iter(distance_vectors.keys()):

        distance_vectors[x] = np.array(distance_vectors[x])
        d_vectors_len += np.alen(distance_vectors[x])

    assert d_vectors_len == length, "Lengths don't match: sequence length = " + str(length) + "; d_vectors length: " + str(d_vectors_len)

    # (distances, frequencies per N words)
    # print distance_vectors, frequencies
    return distance_vectors, durations


def gen_stream(score, sequence, pitch_algorithm, durations_algorithm, score_tempo):

    # keyworded arguments:
    #   window_duration
    #   step
    #   signal_gap

    assert isinstance(pitch_algorithm, PitchAlgorithm) and isinstance(durations_algorithm, DurationsAlgorithm)

    dv, durations = gen_pitch_duration_vectors(sequence, pitch_algorithm, durations_algorithm)

    print len(sequence)

    for x in dv.keys():
        dv[x] = iter(dv[x])

    durations = iter(durations)

    # TODO: integrar esta parte no algoritmo anterior para poupar iteracoes
    scale_len = len(scale.MajorScale().getPitches())
    s = scale.MajorScale()

    part = stream.Part()

    assert isinstance(score_tempo, tempo.MetronomeMark) and score_tempo == score.getElementsByClass(tempo.MetronomeMark)[0]
    # part.insert(0, assigned_instrument)

    print 'Assigning notes and durations from numeric vectors...'

    for l in sequence:

        pitch = dv[l].next()
        d = durations.next()

        if l is not '-':

            n = pitch % scale_len
            n = s.getPitches()[n]
            n = note.Note(n)

            n.duration = duration.Duration(d)

        else:

            n = note.Rest()
            n.duration = duration.Duration(d)

            """else:
                n = pitch % scale_len
                n = s.getPitches()[n]
                n = note.Note(n)
                n.octave = 10

                n.duration = duration.Duration(10.0)"""

        n.addLyric(l)

        assert isinstance(n, Music21Object)
        part.append(n)

    print 'Insering part on score'
    score.insert(0, part)
    print 'Done'


def gen_song(pitch_algorithm, durations_algorithm, dynamics_algorithm, alignment, instruments, piece_length=5000):

    ####### ALIGNMENT HANDLING ##############
    assert (alignment is not None), 'No MSA provided'

    assert alignment \
           and (isinstance(alignment, MultipleSeqAlignment) or
                (isinstance(alignment, str) and os.path.isfile(alignment)))

    if isinstance(alignment, str):

        print 'Reading alignment...'
        aln_file = AlignIO.read(open(alignment, 'rU'), 'clustal')
        aln_file = msa_to_phylip(aln_file)

        print 'Opening phylip file...'
        alignment = AlignIO.read(open(aln_file.split('.')[0] + ".phy"), 'phylip-relaxed')

    assert isinstance(piece_length, int) or isinstance(piece_length, float) and piece_length > 60

    # step = 1 if not pitch_algorithm['n_nucleotides'] else pitch_algorithm['n_nucleotides']

    # piece_length for now is only referring to number of musical elements
    # n_pieces = len(alignment[0]) / (step * piece_length)
    scores = []

    alignment = np.array([[y for y in x] for x in alignment])

    for p in range(0, np.alen(alignment[0]), piece_length):

        score = stream.Score()

        # TODO: Make tempo dynamic
        # score_tempo = tempo.MetronomeMark('Larghissimo', 19)
        score_tempo = tempo.MetronomeMark('adagio', 125)
        score.insert(0, score_tempo)

        print 'Generating pitches and durations...'

        subalignment = alignment[:, p:p+piece_length]

        for i in range(0, np.alen(alignment)):
            # gen_stream(score, subalignment[i], pitch_algorithm, durations_algorithm, instruments[i], score_tempo)
            gen_stream(score, subalignment[i], pitch_algorithm, durations_algorithm, score_tempo)

        print 'Checking if parts have the same total duration...'

        # aligning part durations (score or midi cannot be produced with unequal duration parts)
        for part in score.parts:

            # obtaining highest duration from all parts
            # and aligning with it
            diff = score.highestTime - part.highestTime

            if diff > 0:

                while round(diff, 5) > 0:

                    # minimum = duration.Duration('2048th')
                    n = note.Rest()

                    if diff >= float(0.5):
                        n.duration = duration.Duration(0.5)
                    else:
                        if diff >= MIN_TEMPO:
                            n.duration = duration.Duration(diff)
                        else:
                            n.duration = duration.Duration(MIN_TEMPO)

                    assert n.duration.quarterLength >= MIN_TEMPO

                    part.append(n)
                    diff = score.highestTime - part.highestTime

        # dynamics_vector = gen_dynamics_vector(msa, dynamics_algorithm)
        dynamics_vector = gen_dynamics_vector(subalignment, dynamics_algorithm)

        volumes = dynamics_vector['vol']
        window_size = dynamics_algorithm['window_size']

        score = add_dynamics_to_score(volumes, score, window_size, instruments)
        scores.append(score)

    # parte estatistica e output de ficheiros para @FileWriter
    # retornar score, utilizar dynamics_algorithm, adicionar volumes a score e analisar score
    return scores


def gen_random_seqs(n, MAX, filename):

    short_version = []

    with open('source_sequences/mitochondrion.1.1.genomic.fna', 'rU') as handle:
        seq_iterator = SeqIO.parse(handle, 'fasta')

        i = 0
        for record in seq_iterator:
            if i > MAX:
                break
            else:
                short_version.append(record)
            i += 1

    from random import shuffle
    shuffle(short_version)

    SeqIO.write(short_version[:n], filename, 'fasta')

if __name__ == "__main__":
    # TEST_FILE_1 = SEQ_DIR + '/clustal_11.aln'  # '/mafft_164358.fasta'

    # rnd_file = SEQ_DIR + 'rnd.fasta'
    # gen_random_seqs(8, 40, rnd_file)
    # aln = gen_alignment(rnd_file, n_sequences=8)

    msa = AlignIO.read('output.fasta', 'clustal')
    msa = np.array([np.array(seq) for seq in msa])

    algorithm = DynamicsAlgorithm('entropy', window_size=3000)
    algorithm['criteria'] = 'average'

    d_vector = gen_dynamics_vector(msa, algorithm)

    from math import ceil

    windows = ceil(float(np.alen(msa[0]))/3000)+1
    left = np.arange(windows)
    height = d_vector['entropy']

    print height[-1]
    plt.bar(left, height)
    plt.savefig('fig.png')
    plt.close()

    """
    TEST_FILE_2 = SEQ_DIR + '/mafft_164358.fasta'

    for i in range(10,21):
        gen_song(TEST_FILE_1, duration_mapping=1, block=20, window_duration=i, show_score=False,
                 output_midi='Oryzas_1_' + str(i) + '.mid',audio=True)
        gen_song(TEST_FILE_1, duration_mapping=2, block=20, window_duration=i, show_score=False,
                 output_midi='Oryzas_2_' + str(i) + '.mid',audio=True)

    for i in range(20,26):
        gen_song(TEST_FILE_2, duration_mapping=1, block=20, window_duration=i, show_score=False,
             output_midi='OryzasAmmotheas_1_' + str(i) + '.mid',audio=True)
        gen_song(TEST_FILE_2, duration_mapping=2, block=20, window_duration=i, show_score=False,
             output_midi='OryzasAmmotheas_2_' + str(i) + '.mid',audio=True)"""