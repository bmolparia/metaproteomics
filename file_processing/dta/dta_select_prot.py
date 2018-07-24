#!/usr/bin/env python3

"""
Replace protdb ids in dtaselect file with defline

"""
import argparse
import os
from pymongo import MongoClient


def _run(file, protdb, protdbcoll, host, port, tax, out):
    # if out == True, print to stdout
    # else, write to filename `out`
    client = MongoClient(host, port)
    ProtDB = client[protdb][protdbcoll]

    if args.tax:
        tax_client = MongoClient('wl-cmadmin',27017)
        TaxDB = tax_client.taxonomy.taxonomy
    
    output = []
    with open(args.file) as f:
        line = next(f).rstrip()
        output.append(line)
        while not line.startswith("Locus"):
            line = next(f).rstrip()
            output.append(line)
        for line in f:
            line_split = line.rstrip().split('\t')
            line0 = line_split[0]
            if line0.isnumeric():
                defline = ProtDB.find_one(int(line0))['d']
                d_split = defline.split('|')                
                if args.tax:
                    # This only works with the 081315 version of compil (where the taxid is stored in the defline)
                    try:
                        taxid = d_split[-2]
                        if taxid.startswith("taxon:"):
                            taxid = taxid.replace("taxon:","")
                        taxid = int(taxid)
                        d_split[-2] = TaxDB.find_one({'taxid':taxid})['scientific_name']
                        defline = '|'.join(d_split)
                    except:
                        pass
                if args.ip2:
                    line_split[-1] = d_split[-1]
                    if len(d_split) == 1:
                        line_split[0] = line0 + '||' + defline
                    else:
                        line_split[0] = line0 + '||' + '|'.join(d_split[:-1])
                else:
                     line_split[0] = line0 + '||' + defline
                new_line = '\t'.join(line_split)
                output.append(new_line)
            else:
                output.append(line.rstrip())
    if out == True:
        print("\n".join(output))
        return

    with open(out,'w') as f:
        f.write("\n".join(output))
        return

def run(file, protdb, protdbcoll=None, host='wl-cmadmin', port=27018, tax=False, out=None):
    if not protdbcoll:
        protdbcoll = protdb
    if not out:
        if args.ip2:
            out = os.path.splitext(file)[0]+"_ip2.txt"
        else:
            out = os.path.splitext(file)[0]+"_p.txt"
    _run(file, protdb, protdbcoll, host, port, tax, out)

if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument('--protdb', help='MongoDB ProtDB Database name', type=str, default='ProtDB_072114')
    parser.add_argument('--protdbcoll', help='MongoDB ProtDB Collection name', type=str, default='ProtDB_072114')
    parser.add_argument('--host', help='MongoDB host', type=str, default='wl-cmadmin')
    parser.add_argument('--port', help='MongoDB port', type=int, default=27018)
    parser.add_argument('--tax', help='lookup tax', action = 'store_true')
    parser.add_argument('--ip2', help='format output for IP2 use', action = 'store_true')
    parser.add_argument('-c', help='print to stdout. Otherwise writes to `file`_p.txt', action = 'store_true')
    args = parser.parse_args()
    if not args.c:
        if args.ip2:
            args.c = os.path.splitext(args.file)[0]+"_ip2.txt"
        else:
            args.c = os.path.splitext(args.file)[0]+"_p.txt"

    _run(args.file, args.protdb, args.protdbcoll, args.host, args.port, args.tax, args.c)