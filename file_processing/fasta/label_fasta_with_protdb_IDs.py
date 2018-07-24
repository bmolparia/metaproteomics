#!/usr/bin/env python3

# label_fasta_with_protdb_IDs.py
# 
# simple script to query FASTA deflines against ProtDB and return deflines with ProtDB IDs 
# (assuming each full defline matches exactly one ProtDB document)
#
# usage1: cat my_fasta_file.fasta | python3 label_fasta_with_protdb_IDs.py > my_new_fasta_file.fasta
#
# usage2: cat my_fasta_file.fasta | parallel -j+0 --block 50M --recstart '>' --tmpdir . --pipe python3 label_fasta_with_protdb_IDs.py > my_new_fasta_file.fasta

import sys
from pymongo import MongoClient

def main():
    
    fasta_chunk = sys.stdin

    client = MongoClient('localhost', 27018)
    db = client['ProtDB_072114']
    mongocoll = db['ProtDB_072114']

    for line in fasta_chunk:
        if line[0] == '>':
            print(modify_defline(line, mongocoll), end='')
        else:
            print(line, end='')

def modify_defline(defline, mongocoll):
    real_defline = defline[1:]

    protDB_ID = lookup_protDB_ID(real_defline, mongocoll)

    if protDB_ID:
        modified_defline = ''.join(['>',str(protDB_ID),'||',real_defline])
    else:
        print('!! Couldn\'t lookup ID, exiting', file=sys.stderr)
        sys.exit(1)

    return modified_defline

def lookup_protDB_ID(real_defline, mongocoll):
    
    query = mongocoll.find_one({'d':real_defline.replace('\n','')})

    if query:
        return query['_id']
    else:
        print('!! Error looking up defline: {}'.format(real_defline), file=sys.stderr)
        return None

if __name__ == '__main__':
    main()