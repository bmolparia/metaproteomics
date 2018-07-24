# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:19:08 2016

@author: gstupp
"""

# give a protein cluster a name

from metaproteomics import utils
import os
import re
from itertools import chain
from collections import Counter
import mygene
mg = mygene.MyGeneInfo()
from metaproteomics.analysis import functional_analysis
functionizer = functional_analysis.Functionizer()

from pymongo import MongoClient
interpro = MongoClient(host='wl-cmadmin.scripps.edu').wikidata_src.interpro

"""
If there is any uniprot, use the most common symbol
else:
If there are entrez names:
    query mygene using the most common name
    take the most common symbol
else:
If there are are interpro terms, take the most common one's short name
else:
return the most common data repository name

"""

def name(protein_cluster, verbose=False):
    prot_info = protein_cluster.prot_info
    cluster_prot_ids = protein_cluster.cluster_prot_ids

    if verbose:
        print("=========\n")
        print('\n'.join([x['d'] for x in prot_info]))

    uniprot = [uniprot_fasta_to_prot(x['d']) for x in prot_info if (x['r'].startswith("UniProt") and 'GN=' in x['d'])]
    if uniprot:
        if verbose:
            print("uniprot symbols: {}".format(Counter(uniprot)))
            print("using most common uniprot symbol")
        return Counter(uniprot).most_common(1)[0][0]

    refseq = [refseq_fasta_to_prot(x['d']) for x in prot_info if x['r'] == 'RefSeq']
    if refseq:
        mc_refseq = Counter(refseq).most_common(1)[0][0]
        q = mg.query(mc_refseq, fields='all', species='all')
        if 'hits' in q:
            hits = [x.get('symbol', '') for x in q['hits'] if x.get('symbol', None) != x.get('locus_tag',None)]
            if verbose:
                print("refseq names: {}".format(Counter(refseq)))
                print("mygene symbols: {}".format(Counter(hits)))
                print("using most common refseq name -> gene symbol")
            try:
                return Counter(hits).most_common(1)[0][0]
            except IndexError:
                return ''

    ipa_terms = get_ipa_terms(cluster_prot_ids)
    if ipa_terms:
        ipa_names_lookup = {x['_id']:x['short_name'] for x in interpro.find({'_id':{'$in':ipa_terms}})}
        ipa_names = [ipa_names_lookup.get(x, '') for x in ipa_terms]
        if verbose:
            print("ipa terms: {}".format(Counter(ipa_names)))
        try:
            return Counter(ipa_names).most_common(1)[0][0]
        except IndexError:
            return ''
    r = [x['r'] for x in prot_info]
    if verbose:
        print("repos: {}".format(Counter(r)))
    return Counter(r).most_common(1)[0][0]

def refseq_fasta_to_prot(line):
    last_pipe = line.rindex("|")
    last_bracket = line.rindex("[")
    prot = line[last_pipe+1:last_bracket].strip()
    return prot


def uniprot_fasta_to_prot(line):
    if 'GN=' in line:
        return line[line.index('GN=')+3:].split(' ')[0]
    else:
        return ''

def get_ipa_terms(cluster_prot_ids):
    hashes = functionizer.get_hashes(cluster_prot_ids)
    domain_result = functionizer.get_domains_from_hashes(hashes)
    domain_flat = list(chain(*[list(chain(*[x for x in dr.values()])) for dr in domain_result]))
    ipa_terms = [x['ipa'] for x in domain_flat if 'ipa' in x]
    return ipa_terms


#%%

if __name__ == "__main__":
    BASE = "/home/gstupp/projects/Wolan/cmoon/CM7_CM1E2d56col_unenr123_rawextract/analysis2016_10_17"
    grouped_loci = utils.load(os.path.join(BASE,"grouped_loci_filt_annot_2016_10_17.pkl.gz"))
    for locus in grouped_loci[:20]:
        print(name(locus, verbose=True))
