#!/usr/bin/env python3

# make_fasta_from_ProtDB_IDs_noProtDB.py
# Usage:
# python3 make_fasta_from_ProtDB_IDs_noProtDB.py [ProtDB db name] [ProtDB collection name] --host [mongo host] --port [mongo port]
#
# Looks up all ProtDB IDs (read in from STDIN, one per line in input file)
#
# Can be used as part of a Unix pipeline as follows:
# cat *sqt | grep -P '^M\t1\t' | parallel -j+8 --block 500K --pipe python3 get_protDB_IDs_noProtDBSQT.py "SeqDB_uniprot_human" "SeqDB_uniprot_human" "ProtDB_uniprot_human" "ProtDB_uniprot_human" --host localhost --port 27017 | sort | uniq | python3 make_fasta_from_ProtDB_IDs_noProtDB.py "ProtDB_uniprot_human" "ProtDB_uniprot_human" --host localhost --port 27017 > filtered_DB_firstmatchonly_noProtDB.fasta
#
# Sandip Chatterjee
# 11/1/14

import sys
import argparse
from pymongo import MongoClient

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('db_name', help='MongoDB Database name', type=str)
	parser.add_argument('coll_name', help='MongoDB Collection name', type=str)
	parser.add_argument('--host', help='MongoDB host (mongod or mongos)', type=str)
	parser.add_argument('--port', help='MongoDB port (mongod or mongos)', type=int)
	args = parser.parse_args()

	if args.host:
		host = args.host
	else:
		host = 'localhost'

	if args.port:
		port = args.port
	else:
		port = 27017

	client = MongoClient(host, port)
	db = client[args.db_name]
	coll = db[args.coll_name]

	unique_protDB_IDs = sys.stdin

	for protDB_ID in unique_protDB_IDs:
		#lookup protDB ID, print out FASTA output to STDOUT
		protDB_ID_int = int(protDB_ID.rstrip('\n'))
		query = coll.find_one({'_id':protDB_ID_int})
		# if query: ## if query fails, there's a big problem (ProtDB entry doesn't exist...)
		try:
			print('>',end='')
			if query['d'].startswith('Reverse_'):
				reverse_string = 'Reverse_'
			else:
				reverse_string = ''
			print(reverse_string+str(protDB_ID_int))
			print('\n'.join(split_string_by_n(query['s'], 80)))
		except:
			print('!!!! error retrieving this sequence from ProtDB')
			sys.exit(1)

def split_string_by_n(long_string,n):
	while long_string:
		yield long_string[:n]
		long_string = long_string[n:]

if __name__ == '__main__':
	main()