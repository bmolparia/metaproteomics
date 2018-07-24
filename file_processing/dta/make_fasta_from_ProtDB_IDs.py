#!/usr/bin/env python3

# make_fasta_from_ProtDB_IDs.py
# Usage:
# python3 make_fasta_from_ProtDB_IDs.py [ProtDB db name] [ProtDB collection name] --host [mongo host] --port [mongo port]
#
# Looks up all ProtDB IDs (read in from STDIN, one per line in input file)
#
# For generating a fasta file for use in build_indexDB while preserving numbering

# Given input: 1,80000000
# Should look like:
# >1||sp|P31946|1433B.....
# MTMDKSELVQKAKLAEQAERYDDMAAAMKA....
# >80000000||Reverse_sp|Psg....
# MSDEJFGJUSRRGK

import sys
import argparse
from pymongo import MongoClient

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('db_name', help='MongoDB Database name', type=str)
    parser.add_argument('coll_name', help='MongoDB Collection name', type=str)
    parser.add_argument('--host', help='MongoDB host (mongod or mongos)', type=str, default = 'localhost')
    parser.add_argument('--port', help='MongoDB port (mongod or mongos)', type=int, default = 27017)
    args = parser.parse_args()

    client = MongoClient(args.host, args.port)
    db = client[args.db_name]
    coll = db[args.coll_name]

    unique_protDB_IDs = sys.stdin

    for protDB_ID in unique_protDB_IDs:
        #lookup protDB ID, print out FASTA output to STDOUT
        protDB_ID_int = int(protDB_ID)
        query = coll.find_one({'_id':protDB_ID_int})
        # if query: ## if query fails, there's a big problem (ProtDB entry doesn't exist...)
        try:
            print('>' + str(protDB_ID_int) + '||' + query['d'])
            print(query['s'])
        except:
            print('!!!! error retrieving this sequence from ProtDB')
            sys.exit(1)

if __name__ == '__main__':
    main()