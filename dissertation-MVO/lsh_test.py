import pandas as pd

from random import shuffle
import random
import string

import numpy as np

alphabet = 'acgt'

import datasketch as dk  # MinHash signatures library
from scipy.cluster.hierarchy import linkage, cophenet, fcluster

from music21 import note, stream, duration, tempo
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO


# splits a given sequence from an MSA with alphabet {'A','G','C','T','-'} into shingles
# this technique is based on the k-shingles approach used in document matching algorithms
def split_into_shingles(sequence, k=2):

    length = len(sequence)
    if k < 1 or k > length:
        print 'Invalid parameter k ', k, ' for sequence length ', length
        return None

    # tokenizer
    n_shingles = length - k + 1  # experimental value !!

    shingles = np.zeros((n_shingles, ), dtype="S" + str(k))

    i = 0
    while length - i >= k:

        for s in sequence[i : i + k]:
            shingles[i] += s
        # shingles[i] = [''.join(s) for s in sequence[i : i + k].astype(str)]
        i += 1

    print shingles
    return shingles


def get_minhash_signatures(sets, num_perm=128):
    assert isinstance(sets, np.ndarray) and len(sets.shape) == 3
    hashes = []

    for i in range(0, len(sets)):
        m = dk.MinHash(num_perm=num_perm)

        split_into_shingles()


"""def calc_jaccard_similarities(sets, k=2):
    assert len(set(len(subset) for subset in sets)) == 1

    shingles = np.zeros((len(sets), len(sets[0]) - k + 1, k), dtype="S2")
    minhashes = []

    for i in xrange(len(sets)):

        shingles[i] = split_into_shingles(sets[i], k=k)
        m = dk.MinHash()

        # for s in shingles[i]:
        shingle_str = [''.join(s) for s in shingles[i].astype(str)]
        for s in shingle_str:
            m.update(s.encode('utf-8'))
        minhashes.append(m)

    # print 'Filled shingles', shingles
    print 'K = ', k
    # print len(shingles)

    assert len(sets) == len(minhashes)

    for i in range(1, len(sets)):
        for j in range(0, len(sets)):

            if i != j:
                str1 = [''.join(s) for s in shingles[i].astype(str)]
                str2 = [''.join(s) for s in shingles[j].astype(str)]

                # print 'Shingle:', i, shingles[i], j, shingles[j]

                # print set(str2) & set(str1), len(set(str2) & set(str1))
                # print set(str2) | set(str1), len(set(str2) | set(str1))

                jaccard = minhashes[i].jaccard(minhashes[j])

                print i, j, float(len(set(str2) & set(str1))) / len(set(str2) | set(str1))
                print i, j, jaccard"""


def tokenize_score(score):
    assert isinstance(score, stream.Stream)  # && len(score.parts) <= 1

    assert score.getElementsByClass(tempo.MetronomeMark)

    duration_tokens = np.empty(len(score.getElementsByClass(note.GeneralNote)), dtype="S14")  # dtype=np.dtype()
    note_tokens = np.empty(len(score.getElementsByClass(note.GeneralNote)), dtype="S2")

    i = 0
    for element in score:

        if isinstance(element, note.GeneralNote):
            d = element.seconds
            n = element.name

            duration_tokens[i] = str(d)
            note_tokens[i] = str(n)

            # print duration_tokens[i], note_tokens[i]

            i += 1

    return note_tokens, duration_tokens

""""""""""
# TODO: pensar como organizar estrutura com notas
# 		e como comparar sequencias musicais

p1 = stream.Part()

notes = [note.Note('C'), note.Note('G'), note.Note('C'), note.Note('C#'), note.Note('D'), note.Note('E'), ]
durations = []

import random

for n in notes:
    durations.append(random.uniform(0.15, 0.75))

for d, n in zip(durations, notes):
    n.duration = duration.Duration(d)
    p1.append(n)

notes = [note.Note('C'), note.Note('A'), note.Note('G'), note.Note('C'), note.Note('D'), note.Note('E'), ]

p2 = stream.Part()

for d, n in zip(durations, notes):
    n.duration = duration.Duration(d)
    p2.append(n)

p1.insert(0, tempo.MetronomeMark('adagio', 55))
p2.insert(0, tempo.MetronomeMark('adagio', 55))

p1 = tokenize_score(p1)
p2 = tokenize_score(p2)

print 'Tokens', p1[0], p2[0]"""""""""


def cluster_by_lsh(sets, k=2, num_perm=128):

    # list of 2d ndarrays or 3d ndarray
    assert (isinstance(sets, np.ndarray) and len(sets.shape) == 3) or (isinstance(sets, list) and all(isinstance(x, np.ndarray) and len(x.shape) == 2 for x in sets))                                      # 3d ndarray
    assert isinstance(k, int) and k > 0

    n_pieces = len(sets)
    n_rows = sets.shape[1]
    n_cols = sets.shape[2]

    n_sequence = n_cols + k - 1
    n_shingle_elements = n_sequence - k + 1

    # shingles = np.empty((n_pieces, n_rows * (n_shingle_elements)), dtype="S" + str(k))

    minhashes = []

    # pieces
    shingles_idx = 0
    for p in range(0, n_pieces):

        piece = sets[p]

        minhash = dk.MinHash(num_perm=num_perm)
        shingles = np.empty(n_rows * n_shingle_elements, dtype="S" + str(k))

        # iterating sequences from a region
        for s in range(0, len(piece)):

            # input sequence considering surplus characters
            sequence = np.empty((n_sequence,), dtype="S1")
            sequence[0: n_cols] = piece[s]

            if p != n_pieces - 1: # if we aren't on the last piece:

                next_piece = sets[p+1]
                sequence[n_cols :] = next_piece[s][0 : k-1]  # surplus
            else:

                sequence[n_cols :] = 'Z' # TODO: possivelmente substituir por valor mais provavel

            shingled_sequence = split_into_shingles(sequence, k=k)
            assert len(shingled_sequence) == n_shingle_elements, \
                'Shingled sequence: ' + str(len(shingled_sequence)) + ' and fixed len ' + str(n_shingle_elements)

            # print 'Seq len', len(sequence)
            # print 'Len', len(shingled_sequence), 'Seq', shingled_sequence
            # print 'Shingle len', len(shingles[p][shingles_idx : shingles_idx + n_cols + 1])
            # shingles[p][:, shingles_idx: shingles_idx + len(sequence) - 1] = shingled_sequence
            print shingles_idx + n_cols

            # shingles[p][shingles_idx : shingles_idx + n_shingle_elements] = shingled_sequence
            shingles[shingles_idx: shingles_idx + n_shingle_elements] = shingled_sequence
            shingles_idx += n_shingle_elements

        for word in shingles:
            minhash.update(word)

        minhashes.append(minhash)
        shingles_idx = 0

        # shingle_str = [''.join(s) for s in shingles[piece].astype(str)]
        #for s in shingles[piece]:
        #    minhash.update(s.encode('utf-8'))

    assert len(minhashes) == n_pieces
    print shingles

    distance_matrix = np.empty((n_pieces, n_pieces), dtype=np.float)

    for i in range(0, len(minhashes)):
        for j in range(0, len(minhashes)):

            if i == j:
                distance_matrix[i][j] = 0
            else:

                similarity = minhashes[i].jaccard(minhashes[j])

                if similarity == 0:
                    distance_matrix[i][j] = 1
                else:
                    distance_matrix[i][j] = 1 / similarity

    print distance_matrix

    Z = linkage(distance_matrix)    # todo: ver quais metricas e metodos adequados

    from scipy.cluster.hierarchy import dendrogram

    dendrogram(Z, show_leaf_counts=True)

    import matplotlib.pyplot as plt
    plt.show()

    return fcluster(Z, 0.70)


def calc_jaccard_similarities(sets, k=2, inter_alignments=False):
    assert len(set(len(subset) for subset in sets)) == 1
    assert isinstance(sets, np.ndarray)

    # print sets.shape, len(sets[0])
    # print len(sets[0]) - k + 1
    # print len(sets)

    assert isinstance(inter_alignments, bool)

    minhashes = []

    if inter_alignments:

        assert len(sets.shape) == 3

        shingles = np.empty((sets.shape[0], sets.shape[1] * (sets.shape[2] - k + 1), k), dtype="S2")

        shingle_idx = 0
        set_row_len = sets.shape[2] - k + 1

        for i in range(0, sets.shape[0]):

            m = dk.MinHash()
            for j in range(0, sets.shape[1]):
                # print 'N shingles', len(sets[i][j]) - k + 1

                shingles[i, shingle_idx: shingle_idx + set_row_len] = split_into_shingles(sets[i][j], k=k)
                shingle_idx += set_row_len

            shingle_idx = 0

            shingle_str = [''.join(s) for s in shingles[i].astype(str)]
            for s in shingle_str:
                m.update(s.encode('utf-8'))

            minhashes.append(m)

    # if not inter_alignments:
    else:
        shingles = np.zeros((len(sets), len(sets[0]) - k + 1, k), dtype="S2")

        for i in range(0, len(sets)):

            shingles[i] = split_into_shingles(sets[i], k=k)
            m = dk.MinHash()

            # for s in shingles[i]:
            shingle_str = [''.join(s) for s in shingles[i].astype(str)]
            for s in shingle_str:
                m.update(s.encode('utf-8'))
            minhashes.append(m)

    assert len(sets) == len(minhashes)

    if not inter_alignments:
        n_rows = len(sets) * (len(sets) - 1)
    else:

        n_rows = sets.shape[0] - 1

    # df = pd.DataFrame(data=np.zeros(permutations, 3), index='index', columns=['seq i', 'seq j', 'jaccard'], dtype=np.float)
    # jaccard_dict = {'jaccard' : np.zeros(permutations, dtype=np.float), 'seq i': np.zeros(permutations, dtype=np.int), 'seq j' : np.zeros(permutations, dtype=np.int)}

    jaccard_df = pd.DataFrame(np.empty((n_rows,), dtype=[('i', np.uint8), ('j', np.uint8), ('jaccard', np.float)]))

    row = 0
    for i in range(0, len(sets)):
        for j in range(0, len(sets)):

            if i != j and not ((jaccard_df['i'] == 2) & (jaccard_df['j'] == 5)).any():  # excluding intersections

                str1 = [''.join(s) for s in shingles[i].astype(str)]
                str2 = [''.join(s) for s in shingles[j].astype(str)]

                jaccard = minhashes[i].jaccard(minhashes[j])

                # print i, j, float(len(set(str2) & set(str1))) / len(set(str2) | set(str1))

                jaccard_df['i'][row] = i
                jaccard_df['j'][row] = j
                jaccard_df['jaccard'][row] = jaccard

                row += 1

    # df = pd.DataFrame(data=jaccard_dict)
    return jaccard_df


"""alignment = np.empty((3, 4, 8), dtype="S1")

alignment[0:3] = np.array([[[random.choice(alphabet) for x in y] for y in aln] for aln in alignment])
print alignment"""

import config
msa = AlignIO.read(config.SEQ_DIR + '/clustal3.aln', 'clustal')

msa = np.array([[x for x in y] for y in msa[:, 0:2000]])

msa = np.array([x for x in np.array_split(msa, 5, 1)])

clusters = cluster_by_lsh(msa, k=4)

n_clusters = len(np.unique(clusters))
tempo_range = np.arange(start=35, stop=150, step=(150 - 35) / n_clusters)

print n_clusters
print tempo_range