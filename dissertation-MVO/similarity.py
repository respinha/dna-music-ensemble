"""
    Library to handle MinHash/LSH similarity estimation techniques
"""

import numpy as np
import datasketch as dk  # MinHash signatures library

# import pandas as pd

# alphabet = 'acgt'

"""from music21 import note, stream, duration, tempo
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO"""


# class used to implement textual similarity-based techniques
class SimHandler(object):
    def __init__(self, sets, k=2):
        assert (isinstance(sets, np.ndarray) and len(sets.shape) == 3) \
               or (isinstance(sets, list) and all(isinstance(x, np.ndarray) and len(x.shape) == 2 for x in sets))  # 3d ndarray
        assert isinstance(k, int) and k > 0

        self.k = k
        self.sets = sets

    def cluster_by_similarites(self, threshold=0.7, num_perm=128):
        from scipy.cluster.hierarchy import linkage, cophenet, fcluster

        n_pieces = len(self.sets)
        n_rows = self.sets.shape[1] if isinstance(self.sets, np.ndarray) else len(self.sets[0])

        minhashes = []

        # pieces
        shingles_idx = 0
        for p in range(0, n_pieces):

            piece = self.sets[p]

            minhash = dk.MinHash(num_perm=num_perm)

            n_cols = len(piece[0])
            n_sequence = n_cols + self.k - 1
            n_shingle_elements = n_sequence - self.k + 1

            shingles = np.empty(n_rows * n_shingle_elements, dtype="S" + str(self.k))

            # iterating sequences from a region
            # for s in range(0, len(piece)):
            for s in range(0, n_rows):

                # input sequence considering surplus characters
                sequence = np.empty((n_sequence,), dtype="S1")
                sequence[0 : n_cols] = piece[s]

                if p != n_pieces - 1:  # if we aren't on the last piece:

                    next_piece = self.sets[p + 1]
                    sequence[n_cols:] = next_piece[s][0: self.k - 1]  # surplus
                else:

                    sequence[n_cols:] = 'Z'  # TODO: possivelmente substituir por valor mais provavel

                shingled_sequence = self.__split_into_shingles__(sequence)
                assert len(shingled_sequence) == n_shingle_elements, \
                    'Shingled sequence: ' + str(len(shingled_sequence)) + ' and fixed len ' + str(n_shingle_elements)

                print shingles_idx + n_cols

                shingles[shingles_idx: shingles_idx + n_shingle_elements] = shingled_sequence
                shingles_idx += n_shingle_elements

            for word in shingles:
                minhash.update(word)

            minhashes.append(minhash)
            shingles_idx = 0

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

        Z = linkage(distance_matrix)  # todo: test different metrics

        from scipy.cluster.hierarchy import dendrogram

        dendrogram(Z, show_leaf_counts=True)

        # import matplotlib.pyplot as plt
        # plt.show()
        # plt.savefig('dendrogram_' + str(self.k))

        return fcluster(Z, 0.70)

    def __split_into_shingles__(self, sequence):

        length = len(sequence)
        if self.k < 1 or self.k > length:
            print 'Invalid parameter k ', self.k, ' for sequence length ', length
            return None

        # tokenizer
        n_shingles = length - self.k + 1  # experimental value !!

        shingles = np.zeros((n_shingles,), dtype="S" + str(self.k))

        i = 0
        while length - i >= self.k:

            for s in sequence[i: i + self.k]:
                shingles[i] += s

            i += 1

        print shingles
        return shingles

    def assign_tempos_by_clusters(self, fclusters, tempo_vector):

        # a single tempo value must exist for each
        assert isinstance(tempo_vector, list) or isinstance(tempo_vector, np.ndarray) and len(tempo_vector) >= len(set(fclusters))

        tempo_vector.sort() # sorting if it isn't already sorted
        return np.array([tempo_vector[i] for i in fclusters])