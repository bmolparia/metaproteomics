#!/usr/bin/env python3

"""
Print protdb ids from dtaselect-filter file
"""
import sys
f = sys.stdin
line = next(f).rstrip()
while not line.startswith("Locus"):
    line = next(f).rstrip()
for line in f:
    line_split = line.rstrip().split('\t')
    line0 = line_split[0]
    if line0.isnumeric():
        print(line0)
