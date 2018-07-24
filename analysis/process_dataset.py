"""
Example processing template

For converting dtaselect-filter files into a list of features across all of the samples
A feature is a protein cluster.  Cluster are defined in the clusterdb. 
"""
from metaproteomics import utils
import os
import shelve
from metaproteomics.file_processing import blazmass_tools
from metaproteomics.analysis import build_loci
from pymongo import MongoClient
taxDB = MongoClient('wl-cmadmin', 27017).TaxDB_072114.TaxDB_072114
protDB = MongoClient('wl-cmadmin', 27018).ProtDB_072114.ProtDB_072114
domainDB = MongoClient('wl-cmadmin', 27018).DomainDB_072114.DomainDB_072114
hashDB = MongoClient('wl-cmadmin', 27017).HashDB_072114.HashDB_072114
clusterdb_mongoURI = "mongodb://wl-cmadmin:27017/"
cluster_db_name="071414_ComPIL_forwardonly_0_7"
seqdb_mongoURI="mongodb://wl-cmadmin:27018/"
seqdb='SeqDB_072114'

base = "/mongoc/gstupp/DTASelect/unenriched/UC2UC4UC5"
os.chdir(base)

#%% run on wl-cmadmin
_sample_path = {'UC2_1': "UC2-01/DTASelect-filter_p.txt",
               "UC4_1": "UC4-01/DTASelect-filter_p.txt",
               "UC4_2": "UC4-02/DTASelect-filter_p.txt",
               "UC5_1": "UC5-01/DTASelect-filter_p.txt"}
sample_path = {sample: os.path.join(base,x) for sample,x in _sample_path.items()}
sample_dta = {sample: list(blazmass_tools.dta_select_parser(path, get_tax=True, taxDB=taxDB, protDB=protDB)) for sample,path in sample_path.items()}
with shelve.open("sample_dta.shelve") as f:
    f.update(sample_dta)

#%%
sample_loci = [build_loci.BuildLoci(sample, dta, ppp=2, clusterdb_mongoURI=clusterdb_mongoURI,
                                    cluster_db_name=cluster_db_name, seqdb_mongoURI=seqdb_mongoURI,
                                    seqdb=seqdb) for sample, dta in sample_dta.items()]
utils.save(sample_loci, "sample_loci")
grouped_loci = build_loci.group_across_samples(sample_loci)
utils.save(grouped_loci, "grouped_loci")
build_loci.annotate_function(grouped_loci, protDB=protDB, domainDB=domainDB, hashDB=hashDB)
utils.save(grouped_loci, "grouped_loci_annot")
build_loci.annotate_extra(grouped_loci, protDB=protDB, taxDB=taxDB)
utils.save(grouped_loci, os.path.join(base,"grouped_loci_annot_extra"))