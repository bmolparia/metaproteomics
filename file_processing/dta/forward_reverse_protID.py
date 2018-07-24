"""
If protdbid is forward, print reverse also, and vice verse
Assumes the halfway point between forward and reverse is count() / 2
"""

import sys
import argparse
from pymongo import MongoClient

parser = argparse.ArgumentParser()
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

ProtColl = MongoClient(host, port)[args.protdb_name][args.protdbcoll_name]
half = int(ProtColl.count()/2)

#half = 82817736 # in indexDB / ComPIL
f = sys.stdin
for protID in f:
    protID = int(protID)
    print(protID)
    if protID <= half:
        print(protID + half)
    elif protID > half:
        print(protID - half)
