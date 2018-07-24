# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 15:39:07 2015
build_proteins_from_peptides.py
@author: gstupp

Run DTASelect-filter with p-2. Take all peptides.
Determine all possible protDB IDs

emPAI does not match DTASelect-filter.txt because dtaselect is calculating it: 10^(sequence coverage %) - 1
http://www.mcponline.org/content/4/9/1265.long

"""

from file_processing import blazmass_tools
from itertools import chain
from pymongo import MongoClient
from collections import defaultdict
from numpy import median

def parse(all_peptides, seqDB):
    
    # Get entries in seqDB for each peptide
    # Example seqDB entry: {"_id" : "AAAAAAAAAA", "p" : [ {"r" : "PTG", "d" : true, "i" : 690796, "o" : 475, "l" : "AAK"}] }
    seq_dict_seq = {x['_id']: [i['i'] for i in x['p'] if 'd' not in i] for x in seqDB.find({'_id':{'$in': list(set(iter([x.split('_')[3] for x in all_peptides])))}})}
    seq_dict_seq_rev = {x['_id']: [i['i'] for i in x['p'] if 'd' in i] for x in seqDB.find({'_id':{'$in': list(set(iter([x.split('_')[3] for x in all_peptides])))}})}
    # Make seq_dict match LCStep_Scan_Seq instead of Seq
    seq_dict = {peptide: seq_dict_seq[peptide.split('_')[3]] for peptide in all_peptides}
    seq_dict_rev = {peptide: seq_dict_seq_rev[peptide.split('_')[3]] for peptide in all_peptides}

    # Merge forward & reverse protIDs, make reverse IDs negative
    # Example with both for & rev: Reverse ID: 144540774, peptide: '11_4527_2_RLGVYTR'
    for key in seq_dict.keys(): #seq_dict and seq_dict_rev have the same keys. Just most values are empty in seq_dict_rev
        seq_dict[key].extend([value*-1 for value in seq_dict_rev[key]])

    # Build dict where key is a single protDB_ID, value is set of peptides for that ID
    count = defaultdict(set)
    for peptide, ids in seq_dict.items():
        for id in ids:
            count[id].add(peptide)
    # Frozenset is hashable and can be used as dict keys (below)
    count = dict((key,frozenset(value)) for key, value in count.items())
    
    # Organize into loci
    # key: set of peptides, value: list of protDB_IDs that contain these peptides
    loci = defaultdict(list)
    for (id,p) in count.items():
        loci[p].append(id)
    
    return loci
    
def group_subset(ps):
    ps = sorted(ps, key=lambda x: len(x['psm']))
    new_ps = []
    while len(ps)>0:
        p = ps.pop()
        subs = []
        for x in ps:
            if x['psm'].issubset(p['psm']):
                subs.append(x)                
        ps = [x for x in ps if x not in subs]
        p['subset'] = subs
        new_ps.append(p)
    return new_ps
    
def main(dta_select, ppp = 2, 
         mongo_host = 'wl-cmadmin' , 
         mongo_port = 27018 , 
         seqdb_name = 'SeqDB_072114' ,
         seqdb_coll = 'SeqDB_072114' , 
  		 protdb_name = 'ProtDB_072114' , 
         protdb_coll = 'ProtDB_072114' , 
         coverage = False , 
         emPAI = False):
    
    seqDB = MongoClient(mongo_host, mongo_port)[seqdb_name][seqdb_coll]
    protDB = MongoClient(mongo_host, mongo_port)[protdb_name][protdb_coll]
    
    # Get set of all peptides in the whole file, where a peptide in a particular scan is unique
    # "LCStep_Scan_Seq"
    parser_indexDB = blazmass_tools.dta_select_parser(dta_select, small = False)
    peptides = {'_'.join([x['LCStep'], str(x['Scan']), str(x['ChargeState']), x['AA_Sequence']]): x for x in chain(*[locus['peptides'] for locus in parser_indexDB])}
    
    all_peptides = list(set(peptides.keys()))
        
    loci = parse(all_peptides, seqDB)
    loci = {key: value for (key, value) in loci.items() if len(key) >= ppp} # keep only those with 2 or more peptides
    ps = [{'loci': tuple(value), 'psm': set(key), 'reverse': False, 'Reverse': False} for (key, value) in loci.items()]
    for p in ps:
        if all([x<0 for x in p['loci']]):
            p['reverse'] = True
            p['Reverse'] = True
        p['forward_loci'] = [x for x in p['loci'] if x>0]
        p['reverse_loci'] = [x*-1 for x in p['loci'] if x<0]
        p['loci'] = [abs(x) for x in p['loci']]
        p['peptide_seq'] = set([x.split('_')[3] for x in p['psm']])
        if coverage:
            p['coverage'] = calc_coverage(p, protDB)
        # spectral abundance factor is sum of 'redundancy' counts for all psms in locus,
        # where redundancy is the number of spectra
        p['SAF'] = sum(peptides[x]['Redundancy'] for x in p['psm'])
        
    ps = group_subset(ps)
    
    for p in ps:
        p['parent_forward_loci'] = list(p['forward_loci'])
        p['forward_loci'] = list(chain(*[x['forward_loci'] for x in p['subset']])) + list(p['forward_loci'])
        p['all_loci'] = list(chain(*[x['loci'] for x in p['subset']])) + list(p['loci'])
        if emPAI:
            p['emPAI'] = calc_emPAI(p, seqDB)
    return ps


def calc_coverage(p, protDB):
    coverage = dict()
    sequences = {x['_id']:x['s'] for x in protDB.find({'_id': {'$in': p['loci']}})}
    peptides = p['peptide_seq']    
    for key, s in sequences.items():
        coverage[key] = len(set(list(chain(*[list(range(s.index(pep),s.index(pep) + len(pep))) for pep in peptides])))) / len(s)
    return coverage

def calc_emPAI(p, seqDB):
    numSeq = len(p['peptide_seq'])
    return median([10**(numSeq/seqDB.find({'p.i':locus}).count())-1 for locus in p['parent_forward_loci']])


#[len(x['s']) for x in protDB.find({'_id':{'$in': p['forward_loci']}})]

#dta_select = '/home/gstupp/01_2015_mass_spec/H1_11082014/1108_Gly1_2014_12_15_15_29205/dtaselect_results_sfp0.01_p2/DTASelect-filter.txt'




