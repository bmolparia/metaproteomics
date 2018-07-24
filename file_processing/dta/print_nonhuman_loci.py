#!/usr/bin/env python3
"""
Tool to see what nonhuman loci are found in a sample.
Works on dtaselect-filter file.
Extendable to any taxid, defaults to 9606

Filter out loci with some additional filters: number or peptides, min locus coverage.
Filters out loci that do not have `min peptide` peptides that are NOT in the `target taxid` 
i.e. peptides that never show up in human loci

print in pandas like df output

"""


from metaproteomics.file_processing import blazmass_tools
from metaproteomics.analysis import sample_taxonomy
import os
import argparse
from itertools import chain
import pandas as pd
from collections import Counter
from pymongo import MongoClient
seqDB = MongoClient("wl-cmadmin", 27018).SeqDB_072114.SeqDB_072114

def main():

    parser = argparse.ArgumentParser(description="filter a dta-select filter file")
    parser.add_argument('file', help='path to dta-select filter file', type=str)
    parser.add_argument('-p', '--peptides', help='Min num of pep per prot', type=int, default=2)
    parser.add_argument('-c', '--coverage', help='Min locus coverage', type=float, default=0.1)    
    parser.add_argument('-t', '--taxid', help='target taxid. default: human (9606)', type=int, default=9606)
    parser.add_argument('--brief', help='brief output', action='store_true')
    parser.add_argument('--taxdb-name', default="TaxDB_072114")
    args = parser.parse_args()
    args.file = os.path.expanduser(args.file)
    dta_parser = blazmass_tools.dta_select_parser(args.file)
    st = sample_taxonomy.sample_taxonomy(taxdb_name = args.taxdb_name)
    
    target_loci = []
    non_loci = []
    for locus in dta_parser:
        if locus['reverse']:
            continue
        if len(locus['peptides']) < args.peptides:
            continue
        if not any([x['Sequence Coverage'] >= args.coverage for x in locus['loci']]):
            continue
        taxids = st.get_tax_from_prot(locus['forward_loci'])
        if not taxids:
            continue
        if args.taxid in taxids:
            target_loci.append(locus)
        else:
            non_loci.append(locus)
    
    for locus in non_loci:
        locus_peptides = set([peptide['unmod_peptide'] for peptide in locus['peptides']])
        # Are all of the peptides in this locus also in human?
        # value is # of peptides that could match to a taxID (key)        
        tax_counter = Counter(chain(*[st.get_tax_from_prot([x['i'] for x in seqDB.find_one(peptide)['p']]) for peptide in locus_peptides]))
        if args.taxid in tax_counter and len(locus_peptides) - tax_counter[args.taxid] <= args.peptides:
            # If the number of nonhuman ONLY peptides is less than args.peptides:
            continue
        print(locus['name'])
        if not args.brief:
            df = pd.DataFrame(locus['loci'])
            columns = pd.Series(list(set(df.columns) - set(['Descriptive Name','Validation Status','pI'])))
            print(df[columns].to_string())
            df = pd.DataFrame(locus['peptides'])
            df = df.drop_duplicates("unmod_peptide")
            columns = pd.Series(list(set(df.columns) - set(['mods','aa_sequence','TotalIntensity','isModified','is_modified','lc_step','scan'])))
            print(df[columns].to_string())
            print("Unique peptide count: " + str(len(locus_peptides) - tax_counter[args.taxid]))
        print("------------------------------")
            
if __name__ == "__main__":
    main()


