from Bio import AlignIO
#from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo

from StringIO import StringIO
from Bio import SeqIO

import time
import subprocess
import sys

muscle = r"/usr/local/bin/muscle3.8.31_i86linux64"

muscle_cline = MuscleCommandline(muscle,input=r"source_sequences/short_version3.fna", out=r"source_sequences/muscle.fna", clwstrict=True)

t0 = time.time()
#child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))

stdout, stderr = muscle_cline()
#align = AlignIO.read(stdout, "fasta")
#print align
#stdout, stderr = muscle_cline()
print (time.time() - t0)/60

"""tree = Phylo.read('short_version.dnd', 'newick')
print tree

Phylo.draw_ascii(tree)"""
