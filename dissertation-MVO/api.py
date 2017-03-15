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

from Bio import Cluster

import numpy as np

import os
import sys

from Bio.Phylo.TreeConstruction import DistanceCalculator
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster

from matplotlib import pyplot as plt

from music21 import note, scale, instrument, stream
from music21 import Music21Object

default_mapping = 'default'

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
SEQ_DIR = CURR_DIR + "/source_sequences"

aln = lambda n: SEQ_DIR + '/mafft_' + str(n) + '.fasta'


GLOBALS = {'MEME_URL' : 'http://meme-suite.org/opal2/services/MEME_4.11.2',
           'SEQUENCES' :  SEQ_DIR + '/mitochondrion.1.1.genomic.fna',
           'SUPPORTED ALGORITHMS' : ['clustal', 'mafft', 'muscle'],
           'VALID_MAPPING' : [default_mapping],
           'ALPHABET' : ['A', 'C', 'G', 'T', '-']}


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
def gen_alignment(seq_vector=None, n_sequences=None, algorithm='clustal'):
    assert seq_vector is not None or n_sequences is not None, \
        'Both arguments are None (sequence vector and number of sequences)'

    assert isinstance(seq_vector, list) or isinstance(n_sequences, int), \
        'Either one of two must be provided: sequence vector or number of sequences'

    assert algorithm in GLOBALS['SUPPORTED ALGORITHMS'], \
        'Algorithm does not match any of the currently supported MSA algorithms'

    iterable = SeqIO.parse(open(GLOBALS['SEQUENCES'], 'rU'), 'fasta')

    seq_file = 'pre_alignment.fna'

    if seq_vector is not None:
        sequences = (r for r in iterable if r.description.split('|')[-1] in seq_vector)
    else:
        sequences = generator_from_iterable(iterable, n_sequences)

    SeqIO.write(sequences, seq_file, 'fasta')


    try:
        t0 = time.time()
        if algorithm == 'clustal':
            algorithm = 'clustalw2'
            cline = ClustalwCommandline(algorithm,
                                        infile=seq_file,
                                        outfile='source_sequences/clustal_' + str(n_sequences) + '.aln')
        elif algorithm == 'muscle':
            algorithm = r"/usr/local/bin/muscle3.8.31_i86linux64"
            cline = MuscleCommandline(algorithm, input=seq_file,
                                      out='source_sequences/muscle_' + str(n_sequences) + '.fna',
                                      clwstrict=True)
        elif algorithm == 'mafft':
            algorithm = r"/usr/local/bin/mafft"
            cline = MafftCommandline(algorithm,
                                     input=seq_file,
                                     clustalout=True)

            print cline
            stdout, stderr = cline()

            print 'Elapsed time: ' + str((time.time() - t0) / 60)
            print 'Now writing...\n'

            dst_name = 'source_sequences/mafft_' + str(n_sequences) + '.fasta'
            with open(dst_name, "w") as handle:
                handle.write(stdout)
        else:
            print 'Unknown algorithm\n'

        #stdout, stderr = cline()
        #print 'Elapsed time: ' + str((time.time() - t0) / 60)
        #print 'Now writing...\n'

    except:
        print 'Error aligning with ' + algorithm

    if os.path.isfile(seq_file):
        os.unlink(seq_file)


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
        #tree = Cluster.treecluster(distancematrix=dm)
        #tree.scale()

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
def cluster_alignment(aln_file, depth=1):
    print 'Converting from clustal format to phylip...'

    aln_file = msa_to_phylip(aln_file)

    print 'Opening phylip file...'
    alignment = AlignIO.read(open(aln_file, 'rU'), 'phylip-relaxed')

    print 'Retrieving clusters...'
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
            threshold = i + window
        else:
            threshold = len(sequence)

        subset = sequence[i:threshold]

        # grouping frequencies of A,C,G,T  within a window
        distributions = {l: float(len(filter(lambda y: y == l, subset))) / window for l in alphabet}

        for j in range(i, threshold):

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

            duration = 'eighth'
            keys = durations.keys()

            for k in keys:

                #print str(k) + ', ' + str(dist)
                if k > dist:
                    duration = durations[k]
                    break

            frequencies.append(duration)

    frequencies = np.array(frequencies)

    for x, y in vectors.iteritems():
        diff = length - last_occurrence[x]
        vectors[x].append(diff)

    for x in vectors.keys():
        vectors[x] = np.array(vectors[x])

    # (distances, frequencies per N words)

    print vectors, frequencies
    return vectors, frequencies



def gen_stream(score, sequence):

    dv = distance_vectors(sequence)

    print dv

    for x in dv.keys():
        dv[x] = iter(dv[x])

    scale_len = len(scale.MajorScale().getPitches())
    s = scale.MajorScale()

    part = stream.Part()
    part.insert(0, instrument.Violin())

    print sequence

    for l in sequence:

        if l is not '-':

            n = dv[l].next() % scale_len
            n = s.getPitches()[n]

            n = note.Note(n)

        else:

            n = note.Rest()

        assert isinstance(n, Music21Object)

        part.append(n)

    score.insert(0, part)


def main():

    if len(sys.argv) > 1:
        if sys.argv[1] == 'align':
            gen_alignment(n_sequences=2, algorithm='mafft')

    else:
        alignment = AlignIO.read(open(SEQ_DIR + "/mafft_2.fasta"), 'clustal')
        score = stream.Score()

        for record in alignment:
            gen_stream(score, record.seq[:70])

        score.show(app='/usr/bin/mscore')
        score.write('midi',fp=CURR_DIR + '/alignment.mid')


if __name__ == "__main__":

    #main()
    print numeric_vectors(['A','C','T','G', 'G','T','A','G','T','A','C','C','T','A','G'])
    #print word_frequences(['A','C','T','G', 'G','T','A','G','T','A','C','C','T','A','G'])
