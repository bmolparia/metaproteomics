#!/usr/bin/env python3
"""
make_species_fasta_from_ProtDB_IDs.py

Take a list of protDB Ids from stdin.
look up in protdb, get the species
print fasta for all proteins in these species plus reverse proteins

Forward protDB IDs will stay the same. So they will match up with the original database.
Keeping the same protdbIDs for the reverse means we would have to look them up or know how the forward protIDs are
mapped to the reverse, which we do now, but don't always. So I'm setting the reverse protIDs to 1000000000 + forward ID. 
Assuming we won't ever have a db with a billion proteins in it.

Additionally include all input protdb ids that lack species

cat DTASelect-filter_sfp0_10_p1.txt | python3 ~/metaproteomics/file_processing/dta/print_protDB_ids.py | python3 ~/metaproteomics/file_processing/dta/make_species_fasta_from_ProtDB_IDs.py --host wl-cmadmin --port 27018 ProtDB_072114 ProtDB_072114 > filtered_db_species.fasta

"""
import sys
import argparse

from pymongo import MongoClient


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('db_name', help='MongoDB Database name', type=str)
    parser.add_argument('coll_name', help='MongoDB Collection name', type=str)
    parser.add_argument('--host', help='MongoDB host (mongod or mongos)', type=str, default = 'localhost')
    parser.add_argument('--port', help='MongoDB port (mongod or mongos)', type=int, default = 27018)
    args = parser.parse_args()
    
    client = MongoClient(args.host, args.port)
    db = client[args.db_name][args.coll_name]
    
    ids = list(set(int(x) for x in sys.stdin.read().split()))
    forw_ids = [x['_id'] for x in db.find({"_id":{'$in':ids}}) if not x['d'].startswith("Reverse_")]
    # list of all organisms for these protIDs
    organisms = list(set([x.get('o',None) for x in db.find({"_id":{'$in':forw_ids}})]) - {None})

    # Grab all forward prot IDs of these organisms, then union with all forw ids to include ones with no species
    new_protIDs = list(set(x['_id'] for x in db.find({"o":{'$in':organisms}}) if not x['d'].startswith("Reverse_")) | set(forw_ids))
    
    for chunk in chunks(new_protIDs, 200000): 
        query = db.find({"_id":{'$in':chunk}})
        for doc in query:
            print('>' + str(doc["_id"]) + '||' + doc['d'])
            print(doc['s'])
            print('>' + str(1000000000+doc["_id"]) + '||Reverse_' + doc['d'])
            print(doc['s'][::-1])
