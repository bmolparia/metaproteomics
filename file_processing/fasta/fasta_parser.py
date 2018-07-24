#!/usr/bin/env python3

import sys

def main():
	with open(sys.argv[1]) as f:
		for record in parser(f):
			print(record)

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