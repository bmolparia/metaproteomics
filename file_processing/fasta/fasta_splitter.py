# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 09:18:53 2014

@author: Greg
"""
import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('infile', type = str)
parser.add_argument('size', type = float, help="Size in GB")
args = parser.parse_args()

filename = os.path.splitext(args.infile)[0]
chunk_num = 0
chunk_name = filename + '.fasta.' + str(chunk_num)

chunk_handle = open(chunk_name,'w')
myWriter = SeqIO.FastaIO.FastaWriter(chunk_handle)
myWriter.write_header()

for num,f in enumerate(SeqIO.parse(args.infile, 'fasta')):
    myWriter.write_record(f)
    if chunk_handle.tell()>1000000000*args.size:
        chunk_handle.close()
        chunk_num+=1
        chunk_handle = open(filename + '.fasta.' + str(chunk_num),'w')
        myWriter = SeqIO.FastaIO.FastaWriter(chunk_handle)
        myWriter.write_header()

chunk_handle.close()