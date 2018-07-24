#!/usr/bin/env python3

# remove_reverse_fasta_records.py
#
# reads FASTA file on STDIN, outputs proteins to STDOUT that don't begin with 'Reverse'
#
# made for the ProtDB FASTA defline format:
# >209||gi|430025802|ref|YP_007195275.1| VP2 [Cebus albifrons polyomavirus 1]
# >82818599||Reverse_gi|446730281|ref|YP_007173651.1| phage protein-Restriction nuclease superfamily [Lactobacillus phage phiAQ113]
#
# 
# Sandip Chatterjee

import sys

def main():
	for record in parser(sys.stdin):
		if not record['defline'].split('||')[1].startswith('Reverse_'):
			print(assemble_fasta_record(record['defline'], record['sequence']), end='')

def parser(fasta_file_handle):

	defline, sequence = '', []
	for line in fasta_file_handle:
		if line[0] == '>':
			if defline:
				yield {'defline': defline, 'sequence': ''.join(sequence)}
			defline, sequence = line[1:].rstrip('\n'), []
		else:
			sequence.append(line.rstrip('\n'))

	if defline:
		yield {'defline': defline, 'sequence': ''.join(sequence)}

def split_string_by_n(long_string,n):

	'''
	splits a (long) string into n parts, returning one part at a time.

	(generator function)
	'''

	while long_string:
		yield long_string[:n]
		long_string = long_string[n:]

def assemble_fasta_record(defline, full_sequence, split_num=80):

	'''
	returns a FASTA record as a single string (with newline characters embedded).
	FASTA sequence can be split on to multiple lines using split_num (number of characters to include on each line)
	'''

	split_sequence = '\n'.join(split_string_by_n(full_sequence,split_num))  ## inserts a newline after every (split_num) characters of sequence

	return ''.join(['>', defline, '\n', split_sequence, '\n'])

if __name__ == '__main__':
	main()