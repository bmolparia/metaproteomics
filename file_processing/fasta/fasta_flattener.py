#!/usr/bin/env python3

# FASTA flattener
# Prints 'flattened' FASTA file to STDOUT
# Ignores 'Reverse_' proteins

# 10/17/14
# Sandip Chatterjee

import sys

def main():
	with open(sys.argv[1]) as f:
		for record in parser(f):
			full_defline = protdb_id = record['defline']
			full_defline_split = full_defline.split('||')
			protdb_id, defline_text = full_defline_split[0], '||'.join(full_defline_split[1:])
			if not defline_text.startswith('Reverse_'):
				print(record['sequence']+'\t'+protdb_id)

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

if __name__ == '__main__':
	main()