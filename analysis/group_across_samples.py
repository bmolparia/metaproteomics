# -*- coding: utf-8 -*-
"""
metaproteomics/analysis/group_across_samples.py
Created on Tue Feb 24 22:19:26 2015

@author: gstupp

Take several lists of loci (from all samples in study)
Determine which are in common and a unique ID

Example Input:
Sample: |Peptides|Peptides|
1: | a b c | a b d e |
2: | a b | 
3: | a b c d f | a b d e |

Out:
Locus: Samples
abcdf: 1,2,3 (abc,ab,abcdf)
abde: 1,3 (abde,abcde)

Unique ID for a locus across sample is one of: first hash, first protdb id, string cat of psms
   abcdf  abde
1      3     6
2      2  None
3      5     6

"""

'''
# Example
ps1 = [{'peptide_seq':['a','b','c'], 'emPAI': 3, 'id': '588609b8f6d5a4ff250ca9fa4e6629fb'}, {'peptide_seq':['a','b','d','e'], 'emPAI': 6, 'id': 'f71c7d7664a34f41664f4c7cf24234ef'}]
ps2 = [{'peptide_seq':['a','b'], 'emPAI': 2, 'id': 'abc'}]
ps3 = [{'peptide_seq':['a','b','c','d','f'], 'emPAI': 5, 'id': 'abargergcdf'}, {'peptide_seq':['a','b','d','e'], 'emPAI': 6, 'id': 'f71c7d7664a34f41664f4c7cf24234ef'}]

all_ps = [ps1, ps2, ps3]
'''

from itertools import chain
import pandas as pd

def group(all_ps, field = 'SAF'):
    
    # Determine which field to use for quantification
    # 'psm' : # of psms
    # 'emPAI': empai
    # 'coverage': %seq cov
    
    # Add sample label, psm->set
    for sample, sub_p in enumerate(all_ps):
        for locus in sub_p:
            locus['sample']=sample
            locus['peptide_seq'] = set(locus['peptide_seq'])
    
    data = dict()
    # convert from list of lists of loci to a single list of all loci in all samples
    ps = sorted(chain(*all_ps), key=lambda x: len(x['peptide_seq']))
    new_ps = []
    while len(ps)>0:
        p = ps.pop()
        subs = []
        data[p['id']] = [None]*len(all_ps)
        data[p['id']][p['sample']] = len(p[field]) if field == 'psm' else p[field]
        for x in ps:
            if x['peptide_seq'].issubset(p['peptide_seq']):
                subs.append(x)
                data[p['id']][x['sample']] = len(x[field]) if field == 'psm' else x[field]
        ps = [x for x in ps if x not in subs]
        p['group'] = subs
        new_ps.append(p)
    
    for p in new_ps:
        p['samples'] = [x['sample'] for x in p['group']] + [p['sample']]

    return (new_ps, data)

def group_as_df(all_ps, field = 'SAF'):
    
    # Determine which field to use for quantification
    # 'psm' : # of psms
    # 'emPAI': empai
    # 'coverage': %seq cov
    
    # Add sample label, psm->set
    for sample, sub_p in enumerate(all_ps):
        for locus in sub_p:
            locus['sample']=sample
            locus['peptide_seq'] = set(locus['peptide_seq'])
    
    data = pd.DataFrame(index=range(len(all_ps)))
    # convert from list of lists of loci to a single list of all loci in all samples
    ps = sorted(chain(*all_ps), key=lambda x: len(x['peptide_seq']))
    new_ps = []
    while len(ps)>0:
        p = ps.pop()
        subs = []
        data[p['id']] = None
        data[p['id']][p['sample']] = len(p[field]) if field == 'psm' else p[field]
        for x in ps:
            if x['peptide_seq'].issubset(p['peptide_seq']):
                subs.append(x)
                data[p['id']][x['sample']] = len(x[field]) if field == 'psm' else x[field]
        ps = [x for x in ps if x not in subs]
        p['group'] = subs
        new_ps.append(p)
    
    for p in new_ps:
        p['samples'] = [x['sample'] for x in p['group']] + [p['sample']]

    return (new_ps, data)

def generate_id_PSM(p):
    # For a locus, generate an ID, which is the string of PSMs for that locus
    p['id'] = ';'.join(sorted(p['peptide_seq']))
'''
def generate_id_hash(p):
    # For a locus, generate an ID, which is the alphabetically first md5sum of the list of possible proteins for that set of PSMs
    from pymongo import MongoClient
    client = MongoClient('wl-cmadmin', 27017)
    hashDB = client.HashDB_072114.HashDB_072114
    p['id'] = sorted([x['_id'] for x in hashDB.find({'pID':{'$in': list(p['parent_forward_loci'])}}, {'_id': True })])[0]

def generate_id_protID(p):
    # For a locus, generate an ID, which is the numerically first protDB_ID of the list of possible proteins for that set of PSMs
    p['id'] = sorted(p['parent_forward_loci'])[0]
'''

