#!/usr/bin/env python

#	number_fasta_db.py
#	Sandip Chatterjee
#	v2, June 3, 2014
#	
#	For numbering each record of a multi-protein FASTA file (delimiter: || )
#	Ex:
#	>gi|395443651|ref|YP_006383904.1| cobaltochelatase subunit CobN [Pseudomonas putida ND6]
#	to
#	>1||gi|395443651|ref|YP_006383904.1| cobaltochelatase subunit CobN [Pseudomonas putida ND6]
#	
#	Usage:
#	$ cat protein_DB.fasta | ./number_fasta_db.py > numbered_fasta.fasta

import sys

def main():
	
	DB_file = sys.stdin

	counter = 1
	for line in DB_file:
		if line[0] == '>':
			sys.stdout.write('>'+str(counter)+'||'+line[1:])
			counter += 1
		else:
			sys.stdout.write(line)

if __name__ == '__main__':
	main()