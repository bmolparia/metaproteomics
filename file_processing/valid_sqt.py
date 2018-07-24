"""
quick n dirty check if a sqt file is good
"""
import sys
fname = sys.argv[1]

def is_valid_sqt(fname):
    with open(fname, 'rb') as fh:
        last = fh.readlines()[-1].decode()
    
    if not last.startswith('L'):
        return False
    line_split = last.split()
    if len(line_split) != 4:
        return False
    peptide = line_split[3]
    if not peptide.count('.') == 2:
        return False
    return True

if not is_valid_sqt(fname):
    sys.exit(1)
