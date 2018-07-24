#!/usr/bin/env python3

# remove_protDB_IDs.py
#
# really simple script to remove ProtDB IDs from FASTA deflines...
#
# input: FASTA file with deflines like this:
# >18002261||gi|501123921|ref|WP_012173042.1| cytochrome C oxidase [Azorhizobium caulinodans]
#
# output: FASTA file with deflines like this:
# >gi|501123921|ref|WP_012173042.1| cytochrome C oxidase [Azorhizobium caulinodans]
#
# usage: python3 remove_protDB_IDs.py my_fasta_file.fasta

import sys

def main():
    
    with open(sys.argv[1]) as f, open(sys.argv[1].replace('.fasta', '_noprotdbIDs.fasta'), 'w') as g:
        for line in f:
            if line[0] == '>':
                g.write(strip_protDB_ID(line))
            else:
                g.write(line)

    print('Finished')

def strip_protDB_ID(defline):

    if '||' in defline:
        return '>'+''.join(defline.split('||')[1:])
    else:
        return defline

if __name__ == '__main__':
    main()