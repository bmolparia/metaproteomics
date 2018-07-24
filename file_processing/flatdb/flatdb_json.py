#!/opt/applications/python/2.7.1/gnu/bin/python -u

#	flatdb_json.py
#	Sandip Chatterjee
#	v2, June 13, 2014
#	
#	For generating a JSON (peptide) representation of a SORTED peptide flatfile (sorted by peptide sequence or mass)
#	READS FROM STDIN
#	OUTPUTS TO STDOUT
#	
#	Usage:
#	$ ./flatdb_parse.py seqSorted_flatdb_file.flatdb "seq" | ./flatdb_json.py "seq" > seqSorted_JSON_file.json
#	
#	or
#
#	$ ./flatdb_parse.py massSorted_flatdb_file.flatdb "mass" | ./flatdb_json.py "mass" > massSorted_JSON_file.json

import sys
import json

def main():

	try:
		DB_file = sys.stdin
		flatdb_delimiter = sys.argv[1]
	except:
		print "Requires a DB file on sys.stdin and 'seq' or 'mass' argument"
		sys.exit()

	if flatdb_delimiter == 'mass':
		read_massflatfile(DB_file)
	elif flatdb_delimiter == 'seq':
		read_seqflatfile(DB_file)


def read_massflatfile(DB_file):

	current_chunk = []

	line = DB_file.readline()

	while True:
		defline = ''
		current_chunk = []
		protein_dicts = [] ##	list of dicts (one per peptide), to be unpacked before appending to current_chunk
		if not line:
			break
		if line[0] == '>':
			current_chunk.append(line[1:])		##	remove leading '>'
			while True:
				line = DB_file.readline()
				if not line or line[0] == '>':
					break
				current_chunk.append(line)
		print_massjson_chunk(current_chunk)

def read_seqflatfile(DB_file):

	current_chunk = []

	line = DB_file.readline()

	while True:
		defline = ''
		current_chunk = []
		protein_dicts = [] ##	list of dicts (one per peptide), to be unpacked before appending to current_chunk
		if not line:
			break
		if line[0] == '>':
			current_chunk.append(line[1:])		##	remove leading '>'
			while True:
				line = DB_file.readline()
				if not line or line[0] == '>':
					break
				current_chunk.append(line)
		print_seqjson_chunk(current_chunk)

def print_massjson_chunk(current_chunk):

	current_chunk = [line.rstrip('\n') for line in current_chunk]

	json_dict = {}
	chunk_peptides = []
	json_dict['_id'] = int(current_chunk[0].split()[0])

	for line in current_chunk:
		chunk_peptides.append(line.split()[1])

	json_dict['s'] = list(set(chunk_peptides))

	print json.dumps(json_dict)	##	print json string to STDOUT

def print_seqjson_chunk(current_chunk):

	current_chunk = [line.rstrip('\n') for line in current_chunk]

	json_dict = {}
	parent_proteins = []
	json_dict['_id'] = current_chunk[0].split()[1]

	for line in current_chunk:
		parent_dict = {}
		line_split = line.split()
		parent_dict['l'] = line_split[2]
		parent_dict['r'] = line_split[3]
		parent_dict['o'] = int(line_split[5])
		if 'r' in line_split[4]:
			parent_dict['d'] = True
			parent_dict['i'] = int(line_split[4].lstrip('r'))
		else:
			parent_dict['i'] = int(line_split[4])
		parent_proteins.append(parent_dict)

	json_dict['p'] = parent_proteins

	print json.dumps(json_dict)	##	print json string to STDOUT

if __name__ == '__main__':
	main()