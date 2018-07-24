#!/usr/bin/env python3

# get_protDB_IDs_noProtDBSQT.py
# Looks up parent proteins for peptide matches:
# 	- Takes a series of SQT file 'M' lines as input (on STDIN)
#	- Looks up peptide match in SeqDB
#	- Retrieves parent proteins for peptide from SeqDB
#	- Prints ProtDB IDs for parent proteins to STDOUT (one per line)
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
	parser.add_argument('seqdb_name', help='MongoDB SeqDB Database name', type=str)
	parser.add_argument('seqdbcoll_name', help='MongoDB SeqDB Collection name', type=str)
	parser.add_argument('protdb_name', help='MongoDB ProtDB Database name', type=str)
	parser.add_argument('protdbcoll_name', help='MongoDB ProtDB Collection name', type=str)
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

	SeqDB = client[args.seqdb_name]
	SeqColl = SeqDB[args.seqdbcoll_name]

	ProtDB = client[args.protdb_name]
	ProtColl = ProtDB[args.protdbcoll_name]

	for line in sys.stdin:
		linesplit = line.split('\t')
		peptide = linesplit[-2].split('.')[1]
		seqDB_query = SeqColl.find_one({'_id':peptide})
		if seqDB_query:
			parents = seqDB_query['p']
			for parent in parents:
				protDB_query = ProtColl.find_one({'_id':parent['i']})
				if protDB_query:
					print(protDB_query['_id'])
					#print('>'+protDB_query['d'])
					#print(protDB_query['s'])

if __name__ == '__main__':
	main()
