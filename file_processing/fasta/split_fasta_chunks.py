# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 09:18:53 2014

@author: Greg
"""

fasta_file = r"C:\Users\Greg\Documents\su\1.fasta"
from Bio import SeqIO

handle = open(r"C:\Users\Greg\Documents\su\1_1.fasta","w")
myWriter = SeqIO.FastaIO.FastaWriter(handle)
myWriter.write_header()

for num,f in enumerate(SeqIO.parse(fasta_file, 'fasta')):
    if num%100 == 0:
        myWriter.write_record(f)

handle.close()