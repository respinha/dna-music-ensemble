# needed imports
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from matplotlib import pyplot as plt

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import cophenet

from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans

import numpy as np

np.set_printoptions(precision=5, suppress=True)  # suppress scientific float notation

from music21.instrument import StringInstrument, WoodwindInstrument, BrassInstrument, Percussion

from scipy.stats import entropy, itemfreq


families = {0: StringInstrument, 1: WoodwindInstrument,
            2: BrassInstrument, 3: Percussion}

aln = AlignIO.read(open('source_sequences/clustal_5.aln', 'rU'), 'clustal')

calculator = DistanceCalculator(model='trans')
dm = calculator.get_distance(aln)

X = np.array([row for row in dm])
print X
"""

model = KMeans(n_clusters=3, random_state=0)
model.fit(X)

centroids = model.cluster_centers_
labels = model.labels_

for i in range(0, len(aln[0]), 1500):

    window = 1500 if i+1500 < len(aln[0]) else len(aln[0]) - i

    chunk = np.array([[seq[j] for j in range(i, i+window)]
                        for seq in aln])

    print chunk
    print len(chunk[0])

    distance_matrix = calculator.get_distance()"""

"""for j in range(0, len(set(model.labels_))):

    instr = families[j]

    reduced = [i for i in range(0, len(model.labels_)) if model.labels_[i] == j]

    reduced = np.array([X[i] for i in range(0, len(X)) if i in reduced])
    print reduced

    print instr

    from random import shuffle

    sub = instr.__subclasses__()
    shuffle(sub)
    sub = sub[:len(reduced)]

    print sub
    if len(reduced) > 1:

        Z = linkage(reduced)
        mean = np.mean([x[2] for x in Z])

        max_d = mean

        clusters = fcluster(Z, max_d, criterion='distance')
        print clusters

#plt.title('Hierarchical Clustering')

#dendrogram(
#    Z,
#    leaf_rotation=90.,  # rotates the x axis labels
#    leaf_font_size=8.,  # font size for the x axis labels
#)
#plt.show()"""""