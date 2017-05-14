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

from music21 import note, scale, instrument, stream, duration
from music21.instrument import StringInstrument, WoodwindInstrument, BrassInstrument, PitchedPercussion
from music21 import Music21Object
from music21 import midi


CURR_DIR = os.path.dirname(os.path.abspath(__file__))
SEQ_DIR = CURR_DIR + "/source_sequences"
OUTPUT_FILES = CURR_DIR + '/output_files'
MIN_TEMPO = 0.0625  # 64th

aln = lambda n: SEQ_DIR + '/mafft_' + str(n) + '.fasta'

GLOBALS = {'MEME_URL' : 'http://meme-suite.org/opal2/services/MEME_4.11.2',
           'SUPPORTED ALGORITHMS' : ['clustal', 'mafft', 'muscle'],
           'MAPPINGS' : {0: 'NO_SYNC', 1: 'SYNC_BLOCK_DURATION',2:'SYNC_BLOCL_DURATION_DISCRETE',3:'DURATION_WITH_DISTANCES'},
           'ALPHABET' : ['a', 'c', 'g', 't', '-'],
           'MUSE_SCORE' :   '/usr/bin/mscore',
           'FAMILIES' : {0: StringInstrument, 1: WoodwindInstrument,
                                  2: BrassInstrument, 3: PitchedPercussion},
           'SCORES' :   OUTPUT_FILES + '/scores',
           'MIDI' :   OUTPUT_FILES + '/midi',
           'AUDIO'  :   OUTPUT_FILES + '/audio',
           'HIST_DURATIONS' :   OUTPUT_FILES + '/stats/durations',
            'HIST_NOTES' :   OUTPUT_FILES + '/stats/notes',
           'ALIGNMENT_PARAMS' : ['fasta_file', 'seq_vector', 'n_seq', 'algorithm'],
           }

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
    # assert isinstance(depth, int) and depth > 0
    
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


def dynamics_algorithm(msa, window=500, criteria='local', gap_threshold=0.7, n_levels=5):

        # criteria: local, avg ou median entropy
        assert 'window' is not None and window > 0, 'Empty window for dynamics algorithm'

        aln_len = np.alen(msa[0])

        from math import ceil
        n_windows = ceil(float(aln_len)/window)+1
        gaps_below_thresh = np.zeros(n_windows, dtype=np.bool)

        window_idx = 0

        # first iterating through MSA to identify
        # which windows have a percentage of gaps
        # below the established threshold

        for i in range(0,aln_len, window):

            if i + window > aln_len: window = aln_len - i

            local_is_below_threshold = np.zeros((window,), dtype=np.bool)
            for j in range(i, i+window):

                column = np.array([c for c in msa[:,j]])

                n_gaps = float(np.count_nonzero(column == '-'))

                if n_gaps / len(column) < gap_threshold:
                    local_is_below_threshold[j-i] = True

            if np.alen(np.where(local_is_below_threshold)) < 0.7 * window:
                gaps_below_thresh[window_idx] = True

            window_idx += 1

        # n_windows = np.alen(np.where(gaps_below_thresh)[0])

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

        split_info = np.array_split(np.sort(np.unique(entropies)), n_levels)   # splitting info into classes
        volumes = np.linspace(min_vol, max_vol, num=n_levels)                # vector with all possible volumes

        for i in range(0, int(n_windows)):

            for j in range(0, np.alen(split_info)):
                if entropies[i] == -1:
                    dynamics_vector['vol'][i] = -1
                elif entropies[i] <= split_info[j][-1]:

                    dynamics_vector['vol'][i] = volumes[j]
                    break

        return dynamics_vector, n_windows

# returns a tuple containing:
#   - a word distance vector per word (A,C,G,T)
#   - a label array with the assigned duration of each nucleotide
def numeric_vectors(sequence, block=5, **kwargs):

    # kwargs:
    #   - block_duration
    #   - step
    #   - duration_mapping

    assert sequence is not None

    # step is the number of columns/characters that are mapped
    if 'step' in kwargs:
        step = kwargs['step']
        if not isinstance(step, int):
            print 'Invalid step introduced. Using default value 1.'
            step = 1
    else:
        step = 1

    if 'duration_mapping' in kwargs:
        d_mapping = kwargs['duration_mapping']

        assert isinstance(d_mapping, int) and d_mapping in GLOBALS['MAPPINGS'].keys(),\
                                    "Invalid mapping: " + str(d_mapping)
    else:
        d_mapping = 0

    if d_mapping > 0:
        assert 'block_duration' in kwargs, 'No time was assigned for each block of words'

        block_duration = kwargs['block_duration']
        assert isinstance(block_duration, int) or isinstance(block_duration, float)

        if not isinstance(block_duration, float): block_duration = float(block_duration)

    length = len(sequence)

    # auxiliary structures
    distance_vectors = dict()  # keys are nucleotides; values are np arrays
    last_occurrence = dict()  # aux dict used for word distances

    if d_mapping == 0:
        # fixed size array of labels
        durations = np.chararray((length,), itemsize=10)  # aux vector used for durations
    elif d_mapping == 1 or d_mapping == 2:
        # fixed size array of floats
        durations = np.zeros((length,), dtype=np.float)

    """elif d_mapping == 3:
        # TODO: not efficient; reimplement in next week

        print 'Global itemfrequency'
        global_itemfreq = itemfreq(sequence)

        durations = {x[0]: [np.zeros(int(x[1])), 0] for x in global_itemfreq}"""

    for i in range(0, length, block):

        if i + block <= len(sequence):
            boundary = i + block
        else:
            boundary = len(sequence)

        subset = sequence[i:boundary]

        # Frequency/durations depend on selected mapping
        # grouping frequencies of A,C,G,T  within a block
        # distributions = {l: float(len(filter(lambda y: y == l, subset))) / block for l in alphabet}
        counts = dict(itemfreq(subset))

        if d_mapping == 1:

            counts_sum = 0
            for j in range(i, boundary):
                counts_sum += int(counts[subset[j-i]])

        total_time = 0.0
        for j in range(i, boundary, step):

            # word distances
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

            # word frequencies

            local_count = int(counts[subset[j - i]])
            freq = float(local_count) / len(subset)

            if d_mapping == 0:

                duration_labels = {0.25: 'eighth', 0.5: 'quarter', 0.75: 'half', 1.0: 'whole'}
                local_duration = 'eighth'

                keys = duration_labels.keys()

                for k in keys:

                    if k > freq:
                        local_duration = duration_labels[k]
                        break


            # duration-biased algorithms
            elif d_mapping == 1:

                local_duration = float(block_duration) * float(local_count) / float(counts_sum)
                total_time += local_duration

                assert local_duration > MIN_TEMPO, \
                    'Higher tempo required for each subsequence; too short duration was calculated: ' + str(local_duration)

            # duration biased with discrete atributions
            elif d_mapping == 2:

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
                print 'Invalid mapping introduced: ', str(d_mapping)
                raise NotImplementedError

            durations[j] = local_duration

        # if d_mapping > 0 and d_mapping < 3:
        #     total = np.sum([durations[j] for j in range(i, boundary)])

        if d_mapping == 2:

            if round(float(total_time), 5) != round(float(block_duration), 5):

                ratio = float(block_duration) / float(total_time)
                total_time = 0.0
                for j in range(i, boundary):
                    durations[j] *= ratio

                    assert durations[j] >= MIN_TEMPO,\
                        'Higher tempo required for each subsequence; too short duration was calculated: ' + str(durations[j])

                    total_time += durations[j]

        """elif d_mapping == 3:

            if round(float(total_time), 5) != round(float(block_duration), 5):

                ratio = float(block_duration) / float(total_time)
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
        """
        if d_mapping > 0:
            assert round(float(total_time), 5) == round(float(block_duration), 5), \
                'Durations don\'t  match: ' + str(float(total_time)) + ' ' + str(float(block_duration))

    for x, y in iter(distance_vectors.items()):
        diff = length - last_occurrence[x]
        distance_vectors[x].append(diff)

    for x in iter(distance_vectors.keys()):
        distance_vectors[x] = np.array(distance_vectors[x])

    # (distances, frequencies per N words)
    # print distance_vectors, frequencies
    return distance_vectors, durations


def gen_stream(score, sequence, block=10,
               assigned_instr=instrument.Violin(), duration_mapping=0, **kwargs):

    # keyworded arguments:
    #   block_duration
    #   step
    #   signal_gap

    if 'block_duration' in kwargs:
        block_duration = kwargs['block_duration']
        assert isinstance(block_duration, int) or isinstance(block_duration, float)
    else:
        block_duration = -1

    if 'step' in kwargs:
        step = kwargs['step']
        assert isinstance(step, int)
    else:
        step = 1

    if 'signal_gap' in kwargs:
        signal_gap = kwargs['signal_gap']
        if signal_gap is True:
            raise NotImplementedError

    dv, durations = numeric_vectors(sequence, block_duration=block_duration,
                                    block=block, step=step, duration_mapping=duration_mapping)

    # print sequence
    for x in dv.keys():
        # print len(dv[x])
        dv[x] = iter(dv[x])

    # print durations
    durations = iter(durations)

    scale_len = len(scale.MajorScale().getPitches())
    s = scale.MajorScale()

    part = stream.Part()
    print assigned_instr

    part.insert(0, assigned_instr)

    print 'Assigning notes and durations from numeric vectors...'

    c = 0
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


def gen_song(input_alignment=None, alignment_parameters=None, duration_mapping=0, **kwargs):
    # kwargs:
    #   output_midi=None
    #   output_score=None
    #   show_score=Fals,
    #   audio=None,
    #   block=10
    #   block_duration=3

    assert (input_alignment is not None) ^ (alignment_parameters is not None), \
        'Cannot choose alignment file and sequences to align simultaneously'

    if input_alignment is not None:

        print 'Converting from clustal format to phylip...'
        aln_file = msa_to_phylip(input_alignment)
    else:
        assert isinstance(alignment_parameters, dict)

        if 'fasta_file' not in alignment_parameters.keys():
            print 'Missing input file for alignment'
            sys.exit(0)
        seq_file = alignment_parameters['fasta_file']

        if 'algorithm' not in alignment_parameters.keys():
            print 'Missing alignment algorithm'
        algorithm = alignment_parameters['algorithm']

        seq_vector, n_seq = None, None
        if 'seq_vector' not in alignment_parameters.keys():
            if 'n_seq' not in alignment_parameters.keys():
                print 'Missing explicit vector of sequences or maximum number of sequences'
            else:
                n_seq = alignment_parameters['n_seq']
        else:
            seq_vector = alignment_parameters['seq_vector']

        gen_alignment(seq_file, seq_vector=seq_vector, n_sequences=n_seq, algorithm=algorithm)
        aln_file = AlignIO.read(open(seq_file.split('.')[0] + '.aln', 'rU'), algorithm)
        aln_file = msa_to_phylip(aln_file)

    # TODO: insert clustering
    instruments = [instrument.Violin(), instrument.Flute()]
    score = stream.Score()

    print 'Opening phylip file...'
    alignment = AlignIO.read(open(aln_file.split('.')[0] + ".phy"), 'phylip-relaxed')

    print 'Generating score parts...'

    if 'block_duration' not in kwargs:
        block_duration = 4
    else:
        assert isinstance(kwargs['block_duration'], int) and kwargs['block_duration'] > 0
        block_duration = kwargs['block_duration']

    if 'block' not in kwargs:
        block = 10
    else:
        assert isinstance(kwargs['block'], int) #and kwargs['block'] > block_duration
        block = kwargs['block']

    i = 0
    sequence_names = np.chararray((len(alignment),), itemsize=25)

    for record in alignment:

        sequence_names[i] = record.description
        gen_stream(score, record.seq[:400], block=block,
                   assigned_instr=instruments[i], duration_mapping=duration_mapping,
                   block_duration=block_duration)
        i += 1

    print 'Checking if parts have the same time...'
    for part in score.parts:

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

                assert n.duration.quarterLength > MIN_TEMPO
                part.append(n)
                diff = score.highestTime - part.highestTime

    print 'Generating durations histograms'

    if not os.path.isdir(GLOBALS['HIST_DURATIONS']):
        if not os.path.isdir(OUTPUT_FILES + '/stats'):
            os.makedirs(OUTPUT_FILES + '/stats')
        os.makedirs(GLOBALS['HIST_DURATIONS'])

    if input_alignment is not None:

        hist_dir = input_alignment.split('.')[0].split('/')[-1]

    if 'output_midi' in kwargs:
        output_midi = kwargs['output_midi']
        assert isinstance(output_midi, str)

        durations_dir = GLOBALS['HIST_DURATIONS'] + '/' + output_midi
        notes_dir = GLOBALS['HIST_NOTES'] + '/' + output_midi

        if not os.path.isdir(durations_dir):
            os.makedirs(durations_dir)

        if not os.path.isdir(notes_dir):
            os.makedirs(notes_dir)

    if not os.path.isdir(durations_dir + '/' + hist_dir):
        os.makedirs(durations_dir + '/' + hist_dir)

    if not os.path.isdir(notes_dir + '/' + hist_dir):
        os.makedirs(notes_dir + '/' + hist_dir)

    i = 0
    for part in score.parts:

        durations_notes = np.array([(float(x.duration.quarterLength), x.name) for x in part
                                    if not isinstance(x, instrument.Instrument)],
                                   dtype=[('durations', np.float,), ('notes', 'S5')])

        durations = durations_notes['durations'] # first column
        notes = durations_notes['notes']
        # print notes

        plt.hist(durations, color='b')
        plt.savefig(durations_dir + '/' + hist_dir + '/' + sequence_names[i])
        plt.close()

        from collections import Counter
        # import pandas

        counts = Counter(notes)
        # df = pandas.DataFrame.from_dict(counts, orient='index')
        # print df[0]

        note_names = counts.keys()
        note_values = counts.values()

        width = 1.0
        idx = np.arange(len(note_names))
        plt.bar(idx, note_values, width)
        plt.xticks(idx + width * 0.5, note_names)
        plt.savefig(notes_dir + '/' + hist_dir + '/' + sequence_names[i])
        plt.close()

        i += 1

    print 'Checking output parameters...'
    if 'output_score' in kwargs:

        output_score = kwargs['output_score']
        assert isinstance(output_score, str)

        print 'Writing score...'

        """for part in score.parts:
            for n in part:
                # assert float() >= 0.001953125
                if not isinstance(n, instrument.Instrument):
                    print n.quarterLength"""

        # lily already inserts extension in filename
        output_score = output_score[:-4] if output_score.endswith('.pdf') else output_score

        score.write('lily.pdf', fp=GLOBALS['SCORES'] + '/' + output_score)

        to_unlink = GLOBALS['SCORES'] + '/' + output_score
        if os.path.isfile(to_unlink):
            os.unlink(to_unlink)

    if 'output_midi' in kwargs:

        output_midi = kwargs['output_midi']
        assert isinstance(output_midi, str)

        print 'Writing midi...'
        if not output_midi.endswith('.mid'):
            output_midi += '.mid'

        output_midi = GLOBALS['MIDI'] + '/' + output_midi

        f = midi.translate.streamToMidiFile(score)
        f.open(output_midi, attrib='wb')
        f.write()
        f.close()

    if 'show_score' not in kwargs or (isinstance(kwargs['show_score'], bool) and kwargs['show_score'] is True):

        # by default score is shown
        print 'Displaying score...'
        score.show(app=GLOBALS['MUSE_SCORE'])

    if 'audio' in kwargs:
        # print 'Midi to audio conversion not implemented'

        audio = kwargs['audio']
        assert isinstance(audio, bool)

        if audio:
            audio_name = output_midi.split('.')[0] + '.ogg'

            subprocess.call(["timidity",output_midi,"-Ow","-o",audio_name])
            # raise NotImplemented


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

    # TODO: test and debug
    msa = AlignIO.read(SEQ_DIR + '/mafft_164358.fasta', 'clustal')
    msa = np.array([np.array(seq) for seq in msa])

    d_vector, windows = dynamics_algorithm(msa, window=1500)

    left = np.arange(windows)
    height = d_vector['vol']

    plt.bar(left, height)
    plt.show()
    plt.close()

    """
    TEST_FILE_2 = SEQ_DIR + '/mafft_164358.fasta'

    for i in range(10,21):
        gen_song(TEST_FILE_1, duration_mapping=1, block=20, block_duration=i, show_score=False,
                 output_midi='Oryzas_1_' + str(i) + '.mid',audio=True)
        gen_song(TEST_FILE_1, duration_mapping=2, block=20, block_duration=i, show_score=False,
                 output_midi='Oryzas_2_' + str(i) + '.mid',audio=True)

    for i in range(20,26):
        gen_song(TEST_FILE_2, duration_mapping=1, block=20, block_duration=i, show_score=False,
             output_midi='OryzasAmmotheas_1_' + str(i) + '.mid',audio=True)
        gen_song(TEST_FILE_2, duration_mapping=2, block=20, block_duration=i, show_score=False,
             output_midi='OryzasAmmotheas_2_' + str(i) + '.mid',audio=True)"""