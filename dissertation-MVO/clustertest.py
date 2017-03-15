# needed imports
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from matplotlib import pyplot as plt

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import cophenet


from scipy.spatial.distance import pdist, squareform

import numpy as np

np.set_printoptions(precision=5, suppress=True)  # suppress scientific float notation

aln = AlignIO.read(open('source_sequences/clustal_11.aln', 'rU'), 'clustal')

calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(aln)
X = np.matrix([row for row in distance_matrix])

Z = linkage(X)
print Z

#c, coph = cophenet(Z, squareform(X))
#print c,coph


# calculate full dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')

dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
#plt.savefig('dendogram.png')
plt.show()


max_d = 0.2
clusters = fcluster(Z, max_d, criterion='distance')
print clusters