#!/usr/bin/env python3

"""
Replace blank L lines in sqt file
(if you ran blazmass with no seqdb)

From:
S       6332    6332    2       158     nodea1322       1389.72900390625        88.0    9.999999974752427E-7    18293
M       1       1       1385.71 0.0     1.5189252       6.159424        0       30      -.AAGTGDGS(123.123)RSHSPDVVWD(123.123)ANEGVSGEPFKPAPL(111)INPTQNSIGGIPLITLLSPLLNR.- U
L


To:
S       6332    6332    2       158     nodea1322       1389.72900390625        88.0    9.999999974752427E-7    18293
M       1       1       1385.71 0.0     1.5189252       6.159424        0       30      -.AAGTGDGS(123.123)RSHSPDVVWD(123.123)ANEGVSGEPFKPAPL(111)INPTQNSIGGIPLITLLSPLLNR.- U
L	Reverse_4517113	0	IVK.AAGTGDGSRSHSPDVVWDANEGVSGEPFKPAPLINPTQNSIGGIPLITLLSPLLNR.LPA
L	Reverse_4566997	0	IVK.AAGTGDGSRSHSPDVVWDANEGVSGEPFKPAPLINPTQNSIGGIPLITLLSPLLNR.LPA
L	Reverse_4574401	0	IVK.AAGTGDGSRSHSPDVVWDANEGVSGEPFKPAPLINPTQNSIGGIPLITLLSPLLNR.LPA



"""
import argparse
import re
from pymongo import MongoClient

def get_unmod_peptide(peptide_dots):
    peptide = re.findall('\.(.*)\.', peptide_dots)[0]
    is_modified = True if ')' in peptide else False
    if is_modified:
        unmod_peptide = peptide
        for match in re.finditer('\((.*?)\)', peptide):
            unmod_peptide = unmod_peptide.replace(match.group(), '', 1)
        return unmod_peptide
    return peptide

def lookup_L_lines(line):
    peptide_dots = line.split()[9]
    peptide = get_unmod_peptide(peptide_dots)
    doc = seqDB.find_one({'_id':peptide})
    for p in doc['p']:
        if 'd' in p and p['d']:
            print('L\tReverse_' + str(p['i']) + '\t0\t' + p['l'] + '.' + peptide + '.' + p['r'])
        else:
            print('L\t' + str(p['i']) + '\t0\t' + p['l'] + '.' + peptide + '.' + p['r'])
    
    
if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument('--seqdb', help='MongoDB SeqDB Database', type=str, default='SeqDB_072114')
    parser.add_argument('--uri', help='MongoDB host', type=str, default='mongodb://wl-cmadmin:27018')
    parser.add_argument('--skip', help="Skip L lines that already exist", action = "store_true")
    args = parser.parse_args()

    client = MongoClient(args.uri)
    seqDB = client[args.seqdb][args.seqdb]
    with open(args.file) as f:
        try:
            line = next(f).rstrip()
            while True:
                if line.startswith("M") and not args.skip:
                    print(line)
                    lookup_L_lines(line)
                elif line.startswith("M") and args.skip:
                    print(line)
                    next_line = next(f).rstrip()
                    if len(next_line.strip()) == 1:
                        lookup_L_lines(line)
                    else:
                        while next_line.startswith('L'):
                            print(next_line)
                            next_line = next(f).rstrip()
                        line = next_line
                        continue
                elif line.startswith("L"):
                    pass
                else:
                    print(line)
                line = next(f).rstrip()
        except StopIteration:
            pass


    
    
    
    