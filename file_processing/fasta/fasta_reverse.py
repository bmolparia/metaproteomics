#!/usr/bin/env python3

# fasta_reverse.py
# reverses all FASTA protein sequences in input and prints to screen (STDOUT)
#
# usage: cat FASTA1.fasta FASTA2.fasta ... | python3 fasta_reverse.py [-f] [-p] > Forward_Reverse_concat_fasta_filename.fasta
#
# Sandip Chatterjee
# v3, 7/14/14

import sys
import argparse

def main():

	DB_file = sys.stdin

	aparser = argparse.ArgumentParser(usage='cat FASTA1.fasta FASTA2.fasta ... | python3 fasta_reverse.py [-f] [-p] > Forward_Reverse_concat_fasta_filename.fasta')
	aparser.add_argument('-p', '--protdbids', help='Preserve ProtDB IDs in FASTA deflines (>12345||defline) vs. (>defline). Requires MongoDB.', action='store_true')
	aparser.add_argument('-f', '--writeforwardandreverse', help='Print both forward (input) and reverse proteins to STDOUT (default is reverse proteins only)', action='store_true')
	args = aparser.parse_args()

	mongocoll = None
	if args.protdbids:
		try:
			from pymongo import MongoClient
		except:
			print('This requires a working MongoDB ProtDB instance (and pymongo)', file=sys.stderr)
			sys.exit(1)
		client = MongoClient('localhost', 27018)
		db = client['ProtDB_072114']
		mongocoll = db['ProtDB_072114']

	for record in parser(DB_file):
		reversed_record = make_reversed_record(record, args, mongocoll)
		if args.writeforwardandreverse:
			print_fasta_record(record)
		print_fasta_record(reversed_record)

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

def make_reversed_record(fasta_record, args, mongocoll=None):

	if args.protdbids:
		protDB_ID = fasta_record['defline'].split()[0]
		real_defline = ''.join(fasta_record['defline'].split('||')[1:])
		if mongocoll:
			query = mongocoll.find_one({'d': 'Reverse_'+real_defline})
			if query:
				defline = str(query['_id'])+'||'+query['d']
				reversed_sequence = query['s']
			else:
				print('ProtDB query "{}" failed'.format('Reverse_'+real_defline), file=sys.stderr)
				sys.exit(1)
		else:
			print('This requires a working MongoDB ProtDB instance', file=sys.stderr)
			sys.exit(1)
	else:
		defline = 'Reverse_'+fasta_record['defline']
		reversed_sequence = fasta_record['sequence'][::-1]

	return {'defline':defline, 'sequence':reversed_sequence}

def split_string_by_n(long_string,n):
	while long_string:
		yield long_string[:n]
		long_string = long_string[n:]

def print_fasta_record(fasta_record):

	defline, sequence = fasta_record['defline'], fasta_record['sequence']

	split_sequence = '\n'.join(split_string_by_n(sequence,80))  ## inserts a newline after every 80 characters of sequence
	print('>'+defline)
	print(split_sequence)

if __name__ == '__main__':
	main()