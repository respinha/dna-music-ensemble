# needed imports
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from matplotlib import pyplot as plt

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

from scipy.spatial.distance import pdist

import numpy as np

np.set_printoptions(precision=5, suppress=True)  # suppress scientific float notation

aln = AlignIO.read(open('source_sequences/clustal_7.aln', 'rU'), 'clustal')

calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(aln)
X = np.matrix([row for row in distance_matrix])

Z = linkage(X)

print Z

print '##################################'
Z = linkage(X, 'ward')

print Z
#max_d = 4
#clusters = fcluster(Z, max_d, criterion='distance')
#print clusters

#plt.figure(figsize=(10, 8))
#for i in range(0,X.shape[1]-1):

#    plt.scatter(X[:,i], X[:,i+1], c=clusters, cmap='prism')

#plt.show()

#c, coph = cophenet(Z, pdist(X))
#print c
#print coph

