from Bio import SeqIO


short_version = []

i = 0

f = open('descriptions.txt', 'w')
with open('source_sequences/mitochondrion.1.1.genomic.fna', 'rU') as handle:

	seq_iterator = SeqIO.parse(handle, 'fasta')

	for record in seq_iterator:
            f.write(record.description + '\n')
            i += 1

"""with open('short_version.fna', 'w') as output:
	SeqIO.write(short_version, output, 'fasta')"""

print i
