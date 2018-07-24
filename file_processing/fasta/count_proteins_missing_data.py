#!/usr/bin/python

#	count_proteins_missing_data.py
#	
#	Counts the number of proteins missing sequence data (proteins with the 'x' character in their sequence)
#	2/11/14

import sys

def main():
	try:
		DB_file = sys.argv[1]
	except:
		print "Requires one argument, a protein FASTA file. Exiting."
		sys.exit()

	read_fasta(DB_file)

	print "Finished"


def read_fasta(DB_file):

	complete_count = 0
	incomplete_count = 0

	with open(DB_file,'rb') as f:

		line = f.readline()
		
		while True:
			defline = ''
			sequence_lines = []
			# protein_dicts = [] ##	list of dicts (one per peptide), to be unpacked before appending to current_chunk
			if not line:
				break
			if line[0] == '>':
				defline = line[1:].rstrip('\n')		##	remove leading '>' and trailing newline
				while True:
					line = f.readline()
					if not line or line[0] == '>':
						break
					sequence_lines.append(line.rstrip('\n'))

			full_sequence = ''.join(sequence_lines).replace('\n','')

			if 'X' in full_sequence.upper():
				incomplete_count += 1
			else:
				complete_count += 1

	print "Complete protein records:", complete_count
	print "Incomplete protein records:", incomplete_count

if __name__ == '__main__':
	main()