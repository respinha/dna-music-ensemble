from Bio import SeqIO

import sys

short_version = []

MAX = sys.argv[1]
with open('source_sequences/mitochondrion.1.1.genomic.fna', 'rU') as handle:

	seq_iterator = SeqIO.parse(handle, 'fasta')

	i = 0
	for record in seq_iterator:
		if i > MAX:
			break
		else:
			short_version.append(record)
		i += 1

with open('short_version' + str(MAX) + '.fna', 'w') as output:
	SeqIO.write(short_version, output, 'fasta')
