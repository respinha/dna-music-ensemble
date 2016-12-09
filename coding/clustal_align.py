from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo

from StringIO import StringIO
from Bio import SeqIO

import time
import sys

n = sys.argv[1]
cline = ClustalwCommandline('clustalw2', infile='source_sequences/short_version' + str(n) + '.fna', outfile='source_sequences/clustal' + str(n) + '.aln')

t0 = time.time()
stdout, stderr = cline()
print (time.time() - t0)/60

tree = Phylo.read('source_sequences/short_version' + str(n) + '.dnd', 'newick')
print tree

Phylo.draw_ascii(tree)
