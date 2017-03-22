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

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster

import numpy as np

import os
import sys
import random

from Bio.Phylo.TreeConstruction import DistanceCalculator
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster

from matplotlib import pyplot as plt

from music21 import note, scale, instrument, stream, duration
from music21.instrument import StringInstrument, WoodwindInstrument, BrassInstrument, Percussion
from music21 import Music21Object
from music21 import midi

default_mapping = 'default'

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
SEQ_DIR = CURR_DIR + "/source_sequences"
OUTPUT_FILES = CURR_DIR + '/output_files'
aln = lambda n: SEQ_DIR + '/mafft_' + str(n) + '.fasta'


GLOBALS = {'MEME_URL' : 'http://meme-suite.org/opal2/services/MEME_4.11.2',
           'SUPPORTED ALGORITHMS' : ['clustal', 'mafft', 'muscle'],
           'VALID_MAPPING' : [default_mapping],
           'ALPHABET' : ['a', 'c', 'g', 't', '-'],
           'MUSE_SCORE' :   '/usr/bin/mscore',
           'FAMILIES' : {0: StringInstrument, 1: WoodwindInstrument,
                                  2: BrassInstrument, 3: Percussion},
           'SCORES' :   OUTPUT_FILES + '/scores',
           'MIDI' :   OUTPUT_FILES + '/midi',
           'audio'  :   OUTPUT_FILES + '/audio',
           'ALIGNMENT_PARAMS' : ['fasta_file', 'seq_vector', 'n_seq', 'algorithm']
           }


def find_homologous_regions(alignment_file):
    alignment = AlignIO.read(alignment_file, 'clustal')

    summ = AlignInfo.SummaryInfo(alignment)
    expect_freq = {'A': .25, 'G': .25, 'T': .25, 'C': .25}

    from Bio.Alphabet import IUPAC
    from Bio.SubsMat import FreqTable

    e_freq_table = FreqTable.FreqTable(expect_freq, FreqTable.FREQ, IUPAC.unambiguous_dna)

    i = 0
    while i < len(alignment[0]):
        h = summ.information_content(i, i + 20, e_freq_table=e_freq_table, chars_to_ignore=['N', '-'])
        print h
        i += 20


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
def gen_alignment(input_file, seq_vector=None, n_sequences=None, algorithm='clustal', output_file='output'):

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

    print sequences
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

    except:
        print 'Error aligning with ' + algorithm

    #if os.path.isfile(seq_file):
    #    os.unlink(seq_file)


# retrieves a distance matrix from:
#   a) a multiple sequence alignment
#   b) a file containing a multiple sequence alignment
def get_distance_matrix(msa):

    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(msa)

    return np.array([row for row in distance_matrix])


# clusters all sequences in a MSA
def get_clusters_from_alignment(msa, depth, save_dendrogram=False):
    assert msa is not None and isinstance(msa, MultipleSeqAlignment)
    assert isinstance(depth, int) and depth > 0

    print 'Retrieving distance matrix'
    dm = get_distance_matrix(msa)
    if dm is not None:

        print 'Retrieving cluster tree'

        Z = linkage(dm)

        if save_dendrogram:
            plt.title('Hierarchical Clustering Dendrogram')
            plt.xlabel('sample index')
            plt.ylabel('distance')

            dendrogram(
                Z,
                leaf_rotation=90.,  # rotates the x axis labels
                leaf_font_size=8.,  # font size for the x axis labels
            )
            plt.savefig('dendogram.png')

        max_d = 0.01 # TODO: make this dynamic

        clusters = fcluster(Z, max_d, criterion='distance')
        print clusters

        return clusters

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
    assert os.path.isfile(msa)

    out_file = msa.split('.')[0] + '.phy'
    AlignIO.convert(msa, 'clustal', out_file, 'phylip-relaxed')

    return out_file


# returns a tuple containing:
#   - a word distance vector per word (A,C,G,T)
#   - a label array with the assigned duration of each nucleotide
def numeric_vectors(sequence, window=5, step=1):

    vectors = dict()  # keys are nucleotides; values are np arrays

    alphabet = GLOBALS['ALPHABET']
    durations = {0.25: 'eighth', 0.5: 'quarter', 0.75: 'half', 1.0: 'whole'}

    last_occurrence = dict()
    frequencies = []

    length = len(sequence)
    for i in range(0, length, window):

        if i + window <= len(sequence):
            boundary = i + window
        else:
            boundary = len(sequence)

        subset = sequence[i:boundary]

        # grouping frequencies of A,C,G,T  within a window
        distributions = {l: float(len(filter(lambda y: y == l, subset))) / window for l in alphabet}

        for j in range(i, boundary):

            # distances
            letter = sequence[j]

            if letter not in vectors.keys():

                # vectors[letter] = np.zeros(shape=(length))
                vectors[letter] = []
            else:
                diff = j - last_occurrence[letter]
                vectors[letter].append(diff)

            last_occurrence[letter] = j

            #########

            dist = distributions[subset[j - i]]

            local_duration = 'eighth'
            keys = durations.keys()

            for k in keys:

                if k > dist:

                    local_duration = durations[k]
                    break

            #print str(dist) + ', ' + str(local_duration)
            frequencies.append(local_duration)

    frequencies = np.array(frequencies)

    for x, y in vectors.iteritems():
        diff = length - last_occurrence[x]
        vectors[x].append(diff)

    for x in vectors.keys():
        vectors[x] = np.array(vectors[x])

    # (distances, frequencies per N words)

    #print vectors, frequencies
    return vectors, frequencies


def gen_stream(score, sequence, window=5, step=1, assigned_instr=instrument.Violin(),
               signal_gap=False):

    dv, durations = numeric_vectors(sequence, window=window, step=step)

    for x in dv.keys():

        dv[x] = iter(dv[x])

    durations = iter(durations)

    scale_len = len(scale.MajorScale().getPitches())
    s = scale.MajorScale()

    part = stream.Part()
    part.insert(0, assigned_instr)

    for l in sequence:

        pitch = dv[l].next()
        d = durations.next()
        print d

        if l is not '-':

            n = pitch % scale_len
            n = s.getPitches()[n]
            n = note.Note(n)

            n.duration = duration.Duration(d)
        else:

            if not signal_gap:

                n = note.Rest()
                n.duration = duration.Duration(d)

            else:
                n = pitch % scale_len
                n = s.getPitches()[n]
                n = note.Note(n)
                n.octave = 10

                n.duration = duration.Duration(10.0)

        n.addLyric(l)

        assert isinstance(n, Music21Object)
        part.append(n)

    score.insert(0, part)


def gen_song(input_alignment=None, alignment_parameters=None, output_midi=None, output_score=None, show_score=None, audio=None, step=10):

    assert input_alignment is not None ^ alignment_parameters is not None
    assert output_midi is not None or output_score is not None or show_score is not None

    if input_alignment is not None:
        print 'Converting from clustal format to phylip...'
        aln_file = msa_to_phylip(SEQ_DIR + "/" + input_alignment)
    else:
        assert isinstance(alignment_parameters, dict)

        if 'fasta_file' not in alignment_parameters.keys():
            print 'Missing input file for alignment'
            sys.exit(0)
        seq_file = alignment_parameters['fasta_file']

        if 'algorithm' not in alignment_parameters.keys():
            print 'Missing alignment algorithm'
        algorithm = alignment_parameters

        seq_vector, n_seq = None, None
        if 'seq_vector' not in alignment_parameters.keys():
            if 'n_seq' not in alignment_parameters.keys():
                print 'Missing explicit vector of sequences or maximum number of sequences'
            else:
                n_seq = alignment_parameters['n_seq']
        else:
            seq_vector = alignment_parameters['seq_vector']

        gen_alignment(seq_vector=seq_vector, n_sequences=n_seq, algorithm=algorithm)
        aln_file = input_alignment
        pass

    # TODO: insert clustering
    instruments = [instrument.Violin(), instrument.SopranoSaxophone()]
    score = stream.Score()

    print 'Opening phylip file...'
    alignment = AlignIO.read(open(SEQ_DIR + "/" + aln_file.split('.')[0] + ".phy"), 'phylip-relaxed')

    i = 0
    for record in alignment:
        gen_stream(score, record.seq[:400], window=step, assigned_instr=instruments[i], signal_gap=False)
        i += 1
        # print score.highestTime

    for part in score.parts:

        diff = score.highestTime - part.highestTime

        if diff > 0:
            while score.highestTime - part.highestTime > 0:
                n = note.Rest()
                n.duration = duration.Duration(2.0)

                part.append(n)

    if output_score is not None:

        assert isinstance(output_score, str)

        if not output_score.endswith('.pdf'):
            output_score = output_score + '_' + str(step) + '.pdf'

        score.write('lily.pdf', fp=GLOBALS['SCORES'] + '/')  # + '.pdf')
        os.unlink(CURR_DIR + '/' + output_score + '_' + str(step))

    if output_midi is not None:

        assert isinstance(output_midi, str)

        if not output_midi.endswith('.mid'):
            output_midi = output_midi + '_' + str(step) + '.mid'

        output_midi = GLOBALS['MIDI'] + '/' + output_midi

        f = midi.MidiFile()
        f.open(output_midi, attrib = 'wb')
        f.write()
        f.close()

    if show_score is not None:
        score.show(app=GLOBALS['MUSE_SCORE'])

    if audio is not None:
        print 'Audio online conversion not implemented'


if __name__ == "__main__":

    for window in range(10,50,10):
        gen_song()