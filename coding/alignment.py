from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo

from StringIO import StringIO
from Bio import SeqIO

import time

cline = ClustalwCommandline('clustalw2', infile='short_version.fna')

t0 = time.time()
stdout, stderr = cline()
print time.time() - t0

tree = Phylo.read('short_version.dnd', 'newick')
print tree

Phylo.draw_ascii(tree)
