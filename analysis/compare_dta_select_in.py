# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 14:03:08 2015

@author: gstupp

Compare DTASelect with/without -in and python script build_proteins_from_peptides


"""
from file_processing import blazmass_tools
from file_processing.dta import build_proteins_from_peptides
from itertools import chain

# DTA-select
# /mongoa/DTASelect/10_2014_mass_spec/120314_SC_sampleH1sol_HCD35/indexDB_search_noProtDB
#  --quiet --sfp 0.01 -p 2 
dta_select_indexDB = '/home/gstupp/01_2015_mass_spec/120314_SC_sampleH1sol_HCD35/DTASelect-filter.txt'
parser_indexDB = blazmass_tools.dta_select_parser(dta_select_indexDB)
ps_dta = list(parser_indexDB)
ps_dta_loci = set(chain(*[p['forward_loci'] for p in ps_dta]))

# DTA-select
# /mongoc/gstupp/DTASelect/120314_SC_sampleH1sol_HCD35/indexDB_search_noProtDB/dta_in
#  --quiet --sfp 0.01 -p 2 -in 
dta_select_indexDB_in = '/home/gstupp/01_2015_mass_spec/120314_SC_sampleH1sol_HCD35/dta_in/DTASelect-filter.txt'
parser_indexDB_in = blazmass_tools.dta_select_parser(dta_select_indexDB_in)
ps_dta_in = list(parser_indexDB_in)
ps_dta_in_loci = set(chain(*[p['forward_loci'] for p in ps_dta_in]))

# Python script. Use either -in or not, doesn't matter
ps = build_proteins_from_peptides.main(dta_select_indexDB)
ps_loci = set(chain(*[p['forward_loci'] for p in ps]))
ps_ploci = set(chain(*[p['parent_forward_loci'] for p in ps]))


#########
import taxonomy
from pymongo import MongoClient
t = taxonomy.Taxonomy()
client = MongoClient('wl-cmadmin', 27017)
taxDB = client.taxDB.taxDB

for p in ps:
    tax_result = taxDB.aggregate([{'$match':{'_id':{'$in':list(p['parent_forward_loci'])}}},{'$group': {'_id': None, 'taxID': {'$addToSet': '$taxID'}}}])['result']
    p['taxIDs'] = tax_result[0]['taxID'] if tax_result else []
    tax_result = taxDB.aggregate([{'$match':{'_id':{'$in':list(p['forward_loci'])}}},{'$group': {'_id': None, 'taxID': {'$addToSet': '$taxID'}}}])['result']
    p['taxIDs_in'] = tax_result[0]['taxID'] if tax_result else []
    
    try:
        p['LCA'] = t.LCA(p['taxIDs']) if p['taxIDs'] else []
        p['LCA_in'] = t.LCA(p['taxIDs_in']) if p['taxIDs_in'] else []
    except:
        p['LCA'] = []
        p['LCA_in'] = []
    
len([p for p in ps if set(p['taxIDs']) != set(p['taxIDs_in'])])
lca_diff = [p for p in ps if p['LCA'] != p['LCA_in']]

len([p for p in ps if len(p['taxIDs']) == 0 and len(p['taxIDs_in']) > 0])
len([p for p in ps if len(p['taxIDs']) > 0 and len(p['taxIDs_in']) == 0])