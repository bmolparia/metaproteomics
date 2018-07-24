# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 12:41:39 2015

@author: gstupp

Perform GO term enrichment using hashes

Requires:
sudo apt-get install graphviz libgraphviz-dev
sudo easy_install3 fisher (note: don't waste an hour like I did and try to use pip because it doesn't work for some reason)
git clone https://github.com/tanghaibao/goatools.git
sudo python3 setup.py install
wget http://purl.obolibrary.org/obo/go/go-basic.obo -O ~/go/go-basic.obo


This is definitely not the best way to do this:
Ignoring structure of subsets. Each hash is treated independently, when they are clearly not
A "locus"/a "gene" is a set of peptides and all possible subsets of peptides. Can be represented by the numerically lowest protdb_id

"""

HOST = 'wl-cmadmin'

from pymongo import MongoClient
client = MongoClient(HOST, 27017)
taxDB = client.TaxDB_072114.TaxDB_072114
hashDB = client.HashDB_072114.HashDB_072114
domainDB = client.DomainDB_072114.DomainDB_072114

import os
from file_processing.dta import build_proteins_from_peptides
from itertools import chain

from goatools import GOEnrichmentStudy
from goatools.obo_parser import GODag

def setup_association(forward_loci):
    # Take a list of loci and return a set of GO terms for each hash
    hashes = [x['_id'] for x in hashDB.find({'pID':{'$in': list(forward_loci)}}, {'_id': True })]
    domain_result = domainDB.find({'_id':{'$in': hashes}})
    domain_dict = dict([(x['_id'],x['d']) for x in domain_result]) # dictionary where hashes are keys

    assoc = dict()
    for (hash,ds) in domain_dict.items():
        assoc[hash] = set([x for x in chain(*[d['g'] for d in ds if 'g' in d])])
        
    return assoc

def setup_study_pop(forward_loci):
    # Take a list of loci and return a set of hashes
    return [x['_id'] for x in hashDB.find({'pID':{'$in': list(forward_loci)}}, {'_id': True })]


# One of Ana's Sample
study_indexDB = '/home/gstupp/01_2015_mass_spec/H1_11082014/1108_Gly1_2014_12_15_15_29205/dtaselect_results_sfp0.01_p2/DTASelect-filter.txt'
study_ps = build_proteins_from_peptides.main(study_indexDB)
study_loci = set(chain(*[x['forward_loci'] for x in study_ps]))
study = setup_study_pop(study_loci)

# One of Sandip's microbiome samples
pop_indexDB = '/home/gstupp/01_2015_mass_spec/120314_SC_sampleH1sol_HCD35/DTASelect-filter.txt'
pop_ps = build_proteins_from_peptides.main(pop_indexDB)
pop_loci = set(chain(*[x['forward_loci'] for x in pop_ps]))
pop = setup_study_pop(study_loci and pop_loci)

# set up hash -> GO matching
assoc = setup_association(study_loci and pop_loci)


obo_dag = GODag(obo_file=os.path.expanduser("~/go/go-basic.obo"))

study_sub = study[:1000]
g = GOEnrichmentStudy(pop, assoc, obo_dag, alpha=0.05, study=study, methods=["fdr"])
g.print_summary(min_ratio=None, indent=False, pval=None)


