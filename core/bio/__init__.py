from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import MuscleCommandline

import numpy as np
from music21.instrument import PitchedPercussion

from config import GLOBALS
import os
import sys
import random

import time


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

            max_d = 0.01  # TODO: make this dynamic
            clusters = fcluster(Z, max_d, criterion='distance')

        else:
            print 'Invalid cluster algorithm'
            raise NotImplementedError

        if 'algorithm' in kwargs:

            instruments = []

            if 'instruments_pool' in kwargs['algorithm']:
                pass
            else:

                instruments_pool = PitchedPercussion.__subclasses__()
                random.shuffle(instruments_pool)
                instruments_pool = instruments_pool[:nclusters]

                i = 0
                for cluster in clusters:
                    instruments[i] = instruments_pool[cluster]
                    i += 1

                print('Instruments used', instruments);
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
