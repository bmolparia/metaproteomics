#!/usr/bin/env python

#	flatdb_parse.py
#	Sandip Chatterjee
#	v2, June 13, 2014
#	
#	For parsing and preparing a flatdb file to use with GNU parallel and fasta_peptides_json.py
#	
#	Usage:
#	$ ./flatdb_parse.py sorted_flatdb_file.flatdb "mass"
#
#	or
#
#	$ ./flatdb_parse.py sorted_flatdb_file.flatdb "seq"
#	(prints entire flatdb file to STDOUT)

import sys

def main():
	try:
		DB_file = sys.argv[1]
		flatdb_delimiter = sys.argv[2]
	except:
		print "Requires a peptide flatfile"
		print "Correct usage: $ python flatdb_parse.py sorted_flatdb_file.flatdb"
		sys.exit()

	if flatdb_delimiter == 'mass':
		read_flatfile(DB_file,0)
	elif flatdb_delimiter == 'seq':
		read_flatfile(DB_file,1)

def read_flatfile(DB_file,delim):

	current_chunk = []

	with open(DB_file,'rb') as f:

		line = f.readline()
		while True:
			current_chunk = []
			if not line:
				break

			else:
				current_peptide = line.split()[delim]
			while True:
				current_chunk.append(line)
				line = f.readline()
				if not line or line.split()[delim] != current_peptide:
					break

			print_flatdb_chunk(current_chunk)
	
def print_flatdb_chunk(current_chunk):
	
	sys.stdout.write('>')
	for line in current_chunk:
		sys.stdout.write(line)

if __name__ == '__main__':
	main()