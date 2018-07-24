#!/usr/bin/env python3
"""
Make searchgui fasta file from protDB
Takes one argument: the number of forward sequences

>generic|1|gi 526245011 ref YP_008320337.1  terminase small subunit [Paenibacillus phage phiIBB_Pl23]
>generic|1_REVERSED|gi 526245011 ref YP_008320337.1  terminase small subunit [Paenibacillus phage phiIBB_Pl23]
"""
#%%
import sys
LIMIT = int(sys.argv[1])

from pymongo import MongoClient
protDB = MongoClient("wl-cmadmin",27018).ProtDB_072114.ProtDB_072114
n = 1
for rec in protDB.find():
    if n > LIMIT:
        break
    if rec['d'].startswith("Reverse_"):
        continue
    defline = rec['d']
    defline = defline.replace("|"," ")
    fdefline = ">generic|{}|{}".format(n,defline)
    rdefline = ">generic|{}_REVERSED|{}".format(n,defline)
    print(fdefline)
    print(rec['s'])    
    print(rdefline)
    print(rec['s'][::-1])
    n += 1