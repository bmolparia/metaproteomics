"""
Tools for getting taxonomy information about a sample

Input DTASelect-filter - style dictionary

Functions:
Get count and loci information for superkingdoms:
    archaea bacteria eukaryota viroids viruses
Get count and loci for species
"""
from collections import Counter
from itertools import chain
from pymongo import MongoClient
from . import taxonomy

class sample_taxonomy(object):
    def __init__(self, host='wl-cmadmin', port=27017, taxdb_name = "TaxDB_072114",
                 taxonomy_host = 'wl-cmadmin', taxonomy_port=27017):
        # this is where the NCBI taxonomy db is stored
        self.t = taxonomy.Taxonomy(host=taxonomy_host, port=taxonomy_port)
        # this is where the protDB -> taxID mapping is stored
        self.taxDB = MongoClient(host=host, port=port)[taxdb_name][taxdb_name]
        
    
    def lookup_locus_taxonomy(self, dtaselect_loci):
        """
        Adds a 'tax_id' field to dtaselect_loci
        :param dtaselect_loci: a list of loci dicts parsed from a dtaselect-filter file (small=True)
                fields used: 'forward_loci'. must start with protDB_ID
        :return: None. modifies dtaselect_loci in place
        """
        assert type(dtaselect_loci) == list
        for locus in dtaselect_loci:
            if locus['forward_loci'].count("||"):
                protDB_ids = [int(x.split('||')[0]) for x in locus['forward_loci']]
            else:
                protDB_ids = locus['forward_loci']
            locus['tax_id'] = self.get_tax_from_prot(protDB_ids)
    
    def get_tax_from_prot(self,protDB_ids):
        taxIDs_doc = list(self.taxDB.aggregate(
            [{'$match': {'_id': {'$in': protDB_ids}}}, {'$group': {'_id': None, 'taxid': {'$addToSet': '$taxid'}}}]))
        return taxIDs_doc[0]['taxid'] if taxIDs_doc else []
    
    def get_all_tax_at_rank(self, rank="superkingdom"):
        return {x['scientific_name']: x['taxid'] for x in self.t.taxonomy_coll.find({'rank': rank})}
    
    
    def get_count_at_rank(self, dtaselect_loci, rank="superkingdom"):
        """
        fields used: 'tax_id'
        Adds `rank` field to dtaselect_loci in-place
    
        Returns the number of loci that are annotated at the `rank`
        A locus can be counted more than once if it has proteins that are both human and bacterial for instance
    
        For list of ranks see metaproteomics.taxonomy.Taxonomy
    
        :param dtaselect_loci:
        :param rank
        """
        rank_dict = {x['scientific_name']: x['taxid'] for x in self.t.taxonomy_coll.find({'rank': rank})}
        rank_values = set(rank_dict.values())
        sk_counter = Counter()
    
        for locus in dtaselect_loci:
            taxIDs = locus['tax_id']
            if taxIDs:
                locus_superkingdom = set(chain(*self.t.get_lineages(taxIDs))) & rank_values
                sk_counter.update(locus_superkingdom)
                locus[rank] = locus_superkingdom
            else:
                locus[rank] = set()
    
        sk_counter = dict(sk_counter)
        sk_count = {k: sk_counter.get(v, 0) for k, v in rank_dict.items()}
        return sk_count
