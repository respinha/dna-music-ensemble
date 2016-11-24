from Bio import SeqIO


short_version = []

with open('source_sequences/mitochondrion.1.1.genomic.fna', 'rU') as handle:

	seq_iterator = SeqIO.parse(handle, 'fasta')

	i = 0
	for record in seq_iterator:
		if i > 6:
			break
		else:
			short_version.append(record)
		i += 1

with open('short_version.fna', 'w') as output:
	SeqIO.write(short_version, output, 'fasta')
