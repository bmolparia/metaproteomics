#!/usr/bin/env python3

# usage:
# cat massSorted_renumbered_071414_indexDB_reversed.json | parallel -j+0 --block 100M --tmpdir /mongoa/ --pipe python3 make_massdb_histogram.py > make_massdb_histogram.out

import sys
import json

def main():
    for line in sys.stdin:
        json_obj = json.loads(line)
        mass_bin = bin_mass_value(json_obj['_id'])
        num_peptides = len(json_obj['s'])
        print(str(mass_bin)+'\t'+str(num_peptides))

def bin_mass_value(mass_int_x1000):
    ''' 
    want a mass like 510.123 or 511.123 to be binned into '500'
    and a mass like 550.123 should be binned into '600'
    (not perfect, but close... for example, 549.999 will be binned into '500'...)
    '''
    
    return round(mass_int_x1000/100000)*100

if __name__ == '__main__':
    main()