from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo

from StringIO import StringIO
from Bio import SeqIO

import time

def short(i):
    short_version = []

    MAX = i
    with open('source_sequences/mitochondrion.1.1.genomic.fna', 'rU') as handle:

        seq_iterator = SeqIO.parse(handle, 'fasta')

        i = 0
        for record in seq_iterator:
            if i > MAX:
                    break
            else:
                    short_version.append(record)
            i += 1

        with open('source_sequences/short_version' + str(MAX) + '.fna', 'w') as output:
            SeqIO.write(short_version, output, 'fasta')


def align(n_sequences, algorithm='clustal'):

    sequence_file = 'source_sequences/short_version' + str(n_sequences) + '.fna'

    try:
        t0 = time.time()
        if algorithm == 'clustal':
            algorithm = 'clustalw2'
            cline = ClustalwCommandline(algorithm,
                                        infile=sequence_file,
                                        outfile='source_sequences/clustal_' + str(n_sequences) + '.aln')
        elif algorithm == 'muscle':
            algorithm = r"/usr/local/bin/muscle3.8.31_i86linux64"
            cline = MuscleCommandline(algorithm,
                                    input=sequence_file,
                                    out='source_sequences/muscle_' + str(n_sequences) + '.fna',
                                    clwstrict=True)
        elif algorithm == 'mafft':
            algorithm = r"/usr/local/bin/mafft"
            cline = MafftCommandline(algorithm,
                                    input=r"source_sequences/short_version" + str(n_sequences) + r".fna",
                                    clustalout=True)

            print cline
            stdout, stderr = cline()

            print algorithm + ' with ' + str(i) + 'sequences: ' + str((time.time() - t0) / 60)
            print 'Now writing...\n'

            dst_name = 'source_sequences/mafft_' + str(n_sequences) + '.fasta'
            with open(dst_name, "w") as handle:
                handle.write(stdout)
            return 0
        else:
            print 'Unknown algorithm\n'
            return -1

        stdout, stderr = cline()
        print algorithm + ' with ' + str(i) + 'sequences: ' + str((time.time() - t0)/60)
        print 'Now writing...\n'

    except:
        print 'Error aligning with ' + algorithm + ' in iteration ' + str(n_sequences) + '\n'

    ###########################################

def gen_random_seqs(n, MAX):

    import numpy as np
    short_version = []

    with open('source_sequences/mitochondrion.1.1.genomic.fna', 'rU') as handle:
        seq_iterator = SeqIO.parse(handle, 'fasta')

        i = 0
        for record in seq_iterator:
            if i > MAX:
                break
            else:
                short_version.append(record)
            i += 1

    from random import shuffle
    shuffle(short_version)

    SeqIO.write(short_version[:n], 'test.fasta', 'fasta')

gen_random_seqs(10, 40)

"""""
for i in range(3, 12):

    short(i)
    print 'Done short ' + str(i)    
    align(i, 'clustal')
    align(i, 'muscle')
    align(i, 'mafft')

tree = Phylo.read('short_version.dnd', 'newick')
print tree

Phylo.draw_ascii(tree)"""
