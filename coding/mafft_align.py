from Bio import AlignIO
#from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import Phylo

from StringIO import StringIO
from Bio import SeqIO

import time

mafft = r"/usr/local/bin/mafft"
mafft_cline = MafftCommandline('mafft', input=r"source_sequences/short_version.fna")

t0 = time.time()

stdout, stderr = mafft_cline()

print (time.time() - t0)/60 #time in minutes

dst_name = 'source_sequences/mafft.fasta'
with open(dst_name, "w") as handle:
    handle.write(stdout)
align = AlignIO.read(dst_name, "fasta")

"""tree = Phylo.read('short_version.dnd', 'newick')
print tree

Phylo.draw_ascii(tree)"""
