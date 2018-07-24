"""
Example of building a data matrix of protein clusters across samples from DTASelect-filter files

"""

import os
import shelve
import numpy as np
from tqdm import tqdm
from itertools import chain
from sklearn import preprocessing
from metaproteomics import utils
from metaproteomics.analysis import build_loci
from metaproteomics.analysis.DBInfo import DBInfo

if "__file__" in dir():
    BASE = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test")
else:
    BASE = "/home/gstupp/metaproteomics/analysis/build_loci/test"

db_info = DBInfo("compil_mgm")
metadata = build_loci.read_metadata(os.path.join(BASE,"sample_metadata.csv"))

#%% Parse samples
samples = shelve.open(os.path.join(BASE,"samples.shelve"))
for sample_name, sample_info in tqdm(list(metadata.iteritems())):
    samples[sample_name] = build_loci.Sample(sample_name, sample_info.path, db_info, sample_info)
samples.sync()
#%% Build protein clusters
protein_clusters = shelve.open(os.path.join(BASE,"protein_clusters.shelve"))
for name,sample in samples.items():
    protein_clusters[name] = sample.build_protein_clusters()
protein_clusters.sync()
#%%
grouped_loci = build_loci.group_across_samples(chain.from_iterable(protein_clusters.values()))
for locus in tqdm(grouped_loci):
    locus.annotate()
utils.save(grouped_loci, os.path.join(BASE,"grouped_loci.pkl.gz"))

#%% filtering
grouped_loci = [x for x in grouped_loci if x.passes_thresh(metadata, min_samples_per_group=2, group = "biological")]
utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt.pkl.gz"))
build_loci.to_json(grouped_loci, list(samples.keys()), os.path.join(BASE,"sample.json")) # for use with https://bitbucket.org/stuppie/dtaselect-tables
#%% Normalize
nf = build_loci.yates_normalization(samples)
for locus in grouped_loci:
    locus.normalize(nf)
utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt_norm.pkl.gz"))
#grouped_loci = utils.load(os.path.join(BASE,"grouped_loci_filt_norm.pkl.gz"))
#%%
df = build_loci.to_df(grouped_loci, norm=True)

# Log transform
df = np.log(1+df)
# Scale
X = preprocessing.scale(df)

#%%
samples.close()
protein_clusters.close()