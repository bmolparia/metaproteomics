#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:44:13 2015

@author: gstupp, sandip

python3 write_sqt_L_lines.py [-i input_SQT_file.sqt] [-d mongodb_database name] [-c mongodb_collection_name] [--host mongodb_hostname] [--port mongodb_port]

Read in an SQT file. Pull out peptide sequence from M lines, look up in SeqDB
Write L lines. Don't use L lines from orig sqt file


"""

from pymongo import MongoClient
import os
import sys
import glob
import argparse
from multiprocessing import Pool

def main():
    
    if seqDB.find_one() == None:
        raise ValueError("Couldn't connect to seqDB")

    if not os.path.isdir('L_line_corrected'):
        os.makedirs('L_line_corrected')
    
    if args.sqt_input:
        process_SQT_file(args.sqt_input)
    else:
        sqt_files = glob.glob('*.sqt')
        pool = Pool()
        status = pool.imap_unordered(process_SQT_file, sqt_files)
        pool.close()
        pool.join()

        for result in status:
            if not result:
                print('An error occurred -- resubmit')

def process_SQT_file(sqt_file):

    new_sqt_file = 'L_line_corrected'+'/'+sqt_file

    with open(sqt_file) as f, open(new_sqt_file, 'w') as g:
        for line in f:
            if not line.startswith('L'):
                g.write(line)
            if line.startswith('M'):
                peptide = line.split()[9].split('.')[1]
        
                seqDB_result = seqDB.find_one({'_id':peptide})
                for parent in seqDB_result['p']:
                    l_peptide_r = '.'.join([parent['l'], seqDB_result['_id'], parent['r']])
                    if parent.get('d', False):
                        g.write('L\tReverse_{}\t0\t{}\n'.format(parent['i'], l_peptide_r))
                    else:
                        g.write('L\t{}\t0\t{}\n'.format(parent['i'], l_peptide_r))

    print('New file saved to {}'.format(new_sqt_file))

    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--host', help="mongo host", default = 'localhost')
    parser.add_argument('--port', help="mongo port", type = int, default = 27018)
    parser.add_argument('-d', '--mongo-db', help='seq db name', default = 'SeqDB_072114')
    parser.add_argument('-c', '--mongo-coll', help='seq coll name', default = 'SeqDB_072114')
    parser.add_argument('-i', '--sqt_input', help="SQT input file")
    args = parser.parse_args()
    client = MongoClient(host = args.host, port = args.port)
    seqDB = client[args.mongo_db][args.mongo_coll]

    main()