#!/usr/bin/env python3
"""
for converting compil for use with searchgui
>gi|526245011|ref|YP_008320337.1| terminase small subunit [Paenibacillus phage phiIBB_Pl23]|RefSeq|[]

>generic|1|gi 526245011 ref YP_008320337.1  terminase small subunit [Paenibacillus phage phiIBB_Pl23]
>generic|1_REVERSED|gi 526245011 ref YP_008320337.1  terminase small subunit [Paenibacillus phage phiIBB_Pl23]
"""
#%%
from fasta_parser import parser
import sys
for n,rec in enumerate(parser(sys.stdin),1):
    defline = rec['defline']
    if defline.startswith(">"):
        defline.lstrip(">")
    defline = "|".join(defline.split("|")[:-2])
    defline = defline.replace("|"," ")
    fdefline = ">generic|{}|{}".format(n,defline)
    rdefline = ">generic|{}_REVERSED|{}".format(n,defline)
    print(fdefline)
    print(rec['sequence'])    
    print(rdefline)
    print(rec['sequence'][::-1])