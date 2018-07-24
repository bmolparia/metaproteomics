#!/usr/bin/env python3

# usage:
# cat seqSorted_renumbered_071414_indexDB_reversed.json | parallel -j+0 --block 100M --tmpdir /mongoa/ --pipe python3 lookup_reverse_peptides.py > lookup_peptides.out

import sys
import json

def main():
    for line in sys.stdin:
        json_obj = json.loads(line.rstrip('\n'))
        if check_parents_for_decoy(json_obj['p']) == 2:
            print(len(json_obj['_id']))

def check_parents_for_decoy(parents_of_peptide):
    # returns 0 if none of the parents are decoy proteins
    # returns 1 if all of the parents are decoy proteins
    # returns 2 if some (but not all) of the parents are decoy proteins
    
    number_of_decoy_parents = len([parent for parent in parents_of_peptide if 'd' in parent])
    
    if number_of_decoy_parents == 0:
        # no decoy parents
        return 0
    elif number_of_decoy_parents == len(parents_of_peptide):
        # all parents are decoy
        return 1
    else:
        # some (but not all) parents are decoy
        return 2


if __name__ == '__main__':
    main()