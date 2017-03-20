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

aln = AlignIO.read(open('source_sequences/clustal_11.aln', 'rU'), 'clustal')

calculator = DistanceCalculator(model='trans')
distance_matrix = calculator.get_distance(aln)
X = np.array([row for row in distance_matrix])

model = KMeans(n_clusters=3, random_state=0)
model.fit(X)

print model.labels_

centroids = model.cluster_centers_

reduced = [i for i in range(0, len(model.labels_)) if model.labels_[i] == 0]
reduced = np.array([X[i] for i in range(0, len(X)) if i in reduced])

Z = linkage(reduced)
mean = np.mean([x[2] for x in Z])

max_d = mean

clusters = fcluster(Z, max_d, criterion='distance')
print clusters

plt.title('Hierarchical Clustering')

dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.show()

"""Z = linkage(X)
print Z
import pylab

print mean
print type(mean)"""

#plt.figure(figsize=(25, 10))
#plt.title('Hierarchical Clustering Dendrogram')
#plt.title('K-Means Clustering')

#plt.savefig('dendogram.png')
#plt.show()