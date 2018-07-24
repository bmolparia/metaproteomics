#!/usr/bin/env python3
"""
Tool for command line access of metaproteomics.analysis scripts

Looks for files "metadata.csv" and "analysis_params.json" in CWD

# /home/gstupp/projects/Wolan/cmoon/2015_12_01_CM7_CM1E2d56cec_FP_rawextract

"""


import os
import shelve
import numpy as np
import json
from tqdm import tqdm
import argparse
from itertools import chain
from textwrap import dedent
from collections import defaultdict
from sklearn import preprocessing
from metaproteomics import utils
from metaproteomics.analysis import build_loci
from metaproteomics.analysis.DBInfo import DBInfo

BASE = os.getcwd()
params = json.load(open(os.path.join(BASE,"analysis_params.json")))
db_info = DBInfo(params['db_info'])
#%%
def parse_samples(metadata, quiet=False):
    samples = shelve.open(os.path.join(BASE,"samples.shelve"))
    already_done_samples = samples.keys()
    for sample_name, sample_info in tqdm(list(metadata.iteritems()), disable=quiet):
        if sample_name in already_done_samples:
            print("Sample {} already found. Skipping".format(sample_name))
            continue
        samples[sample_name] = build_loci.Sample(sample_name, sample_info.path, db_info, sample_info)
    samples.sync()
    return samples

def build_protein_clusters(samples, quiet=False):
    protein_clusters = shelve.open(os.path.join(BASE,"protein_clusters.shelve"))
    already_done_pc = protein_clusters.keys()
    for sample_name,sample in tqdm(samples.items()):
        if sample_name in already_done_pc:
            print("Protein cluster {} already found. Skipping".format(sample_name))
            continue
        protein_clusters[sample_name] = sample.build_protein_clusters()
    protein_clusters.sync()
    return protein_clusters


def group_across_samples(protein_clusters, sample_pep_quant):
    grouped_loci = build_loci.group_across_samples(chain.from_iterable(protein_clusters.values()), sample_pep_quant, db_info=db_info)
    utils.save(grouped_loci, os.path.join(BASE,"grouped_loci.pkl.gz"), force=True)
    return grouped_loci


def do_filtering(grouped_loci, metadata, f):
    print("Before filtering: {} protein clusters".format(len(grouped_loci)))
    grouped_loci = [x for x in grouped_loci if x.passes_thresh(metadata, min_quant = f['min_quant'], min_samples = f['min_samples'], 
                                           min_samples_per_group = f['min_samples_per_group'], group = f['group'])]
    print("After filtering: {} protein clusters".format(len(grouped_loci)))
    utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt.pkl.gz"), force=True)
    return grouped_loci


def do_annotations(grouped_loci, quiet=False):
    print("Annotating protein clusters")
    build_loci.annotate(grouped_loci, db_info)
    utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt_annot.pkl.gz"), force=True)


def do_normalization(grouped_loci, samples, normalization_type):
    print("Doing normalization")
    norm_func = build_loci.__getattribute__(normalization_type)
    nf = norm_func(samples)
    print("Normalization factors: {}".format(nf))
    for locus in grouped_loci:
        locus.normalize(nf)
    utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt_annot_norm.pkl.gz"), force=True)
        
        
def make_datatables_json(grouped_loci, sample_names, datatables_filename, do_normalization):
    print("Making datatables json")
    build_loci.to_json(grouped_loci, sample_names, os.path.join(BASE, datatables_filename), functionizer=grouped_loci[0].functionizer, norm=do_normalization)


def ssh_datatables_json_wl(datatables_filename):
    #Attempt to ssh this to wl-cmadmin
    import os
    import paramiko
    import getpass
    print("SSHing datatables json to wl-cmadmin")
    k = paramiko.RSAKey.from_private_key_file(os.path.expanduser(os.path.join("~",".ssh","id_rsa")))
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect("wl-cmadmin.scripps.edu", username=getpass.getuser(), pkey=k)
    sftp = c.open_sftp()
    sftp.put(os.path.join(BASE, datatables_filename), os.path.join("/mongoc/dtaselect-tables/static", datatables_filename))
    sftp.close()
    c.close()
    print("http://wl-cmadmin:8000/clustertable/" + os.path.splitext(datatables_filename)[0])
    

def make_df(grouped_loci, use_normalization):
    #Save out pandas Df
    df = build_loci.to_df(grouped_loci, norm=use_normalization)
    df.fillna(0, inplace=True)
    df = np.log(df+1)
    df.to_csv(os.path.join(BASE,"X.csv"))

def do_subtraction(grouped_loci, metadata, params):
    if params['filtering']['do_filtering']:
        thresh = params['filtering']['min_quant']
    
    groupby = params['subtraction']['groupby']
    if groupby:
        # subtract each pair of samples from each other
        for sampletype, this_group in metadata.T.groupby(groupby):
            control = this_group[this_group["control"]].index[0]
            probe = this_group[this_group["control"]==False].index[0]
            for locus in grouped_loci:
                locus.quantification[probe] = locus.quantification.get(probe,0) - locus.quantification.get(control,0)
                if control in locus.quantification:
                    del locus.quantification[control]
    else:
        # subtract all protein clusters in any control sample from all non-control samples
        control_samples = list(metadata.loc[:,metadata.loc[params['subtraction']['control']]==True].columns)
        new_pc = set([x.cluster_id for x in grouped_loci if sum([v for k,v in x.quantification.items() if k in control_samples])<thresh])
        print("#Probe: {}\t#ctrl: {}".format(len(new_pc),len(grouped_loci)-len(new_pc)))
        grouped_loci = [x for x in grouped_loci if x.cluster_id in new_pc]
        
    
    return grouped_loci

def run_pipeline(metadata, params, quiet = False):
        
    samples = parse_samples(metadata, quiet=quiet)
    sample_pep_quant = {sample.sample_name:sample.pep_quant for sample in samples.values()}
    protein_clusters = build_protein_clusters(samples, quiet=quiet)
    grouped_loci = group_across_samples(protein_clusters, sample_pep_quant)

    if params['filtering']['do_filtering']:
        grouped_loci = do_filtering(grouped_loci, metadata, params['filtering'])
    else:
        print("No filtering")
        
    if params['normalization']['do_normalization']:
        do_normalization(grouped_loci, samples, params['normalization']['type'])
    else:
        print("No normalization")
    
    if params['subtraction']['do_subtraction']:
        grouped_loci = do_subtraction(grouped_loci, metadata, params)
        utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt_norm_sub.pkl.gz"), force=True)
    else:
        print("No probe subtraction")
    
    do_annotations(grouped_loci, quiet=quiet)
    sample_names = set(chain(*[locus.quantification.keys() for locus in grouped_loci]))
    make_datatables_json(grouped_loci, sample_names, params['datatables_filename'], params['normalization']['do_normalization'])
    ssh_datatables_json_wl(params['datatables_filename'])
    make_df(grouped_loci, params['normalization']['do_normalization'])


def make_analysis_job(dir_path, queue):
    """make a job file that just runs analyze_dtaselect --run-locally"""
    job_boilerplate = '\n'.join(['#!/bin/bash',
                                 '#PBS -q {}'.format(queue),
                                 '#PBS -l nodes=1:ppn=2',
                                 '#PBS -l walltime=4:00:00',
                                 '#PBS -j oe',
                                 '#PBS -l mem=4gb',
                                 '#PBS -N "{}"'.format(dir_path),
                                 '#PBS -o "{}"'.format(dir_path+'/analysis.$PBS_JOBID')])
    base_job_file = dedent("""
                    echo "################################################################################"
                    echo "Folder: {dir_path}"
                    echo "Running on node: `hostname`"
                    echo "################################################################################"
                    module load python/3.5.1
                    PYTHONPATH=~/lib:$PYTHONPATH
                    source /gpfs/home/gstupp/metaproteomics/venv3.5/bin/activate
                    cd {dir_path}
                    /gpfs/home/gstupp/metaproteomics/scripts/analyze_dtaselect.py --run-locally --quiet
                    """).format(dir_path=dir_path)
    
    job_file_path = os.path.join(dir_path,"analysis.job")           
    with open(job_file_path, 'w') as f:
        print('Writing job file: ' + job_file_path)
        f.write(job_boilerplate + '\n' + base_job_file)
    
    os.system("qsub {}".format(job_file_path))
    

def do_venn(args):
    import os
    import shelve
    from metaproteomics.file_processing import blazmass_tools
    from metaproteomics.analysis import venn
    
    samples = args.samples.split(",")
    venn_type = args.type
    
    # check the samples requested are present
    sample_dta = shelve.open(os.path.join(BASE,"samples.shelve"))
    assert len(sample_dta.keys() & set(samples)) == len(set(samples)), "Samples {} not found".format(set(samples) - sample_dta.keys())    
    
    if venn_type == "peptide":
        sample_peptides = {sample_name:blazmass_tools.get_unmod_peptides(sample.dta_select) for sample_name, sample in sample_dta.items() if sample_name in samples}
        print(sample_peptides.keys())
        venn.make_vens(sample_peptides, samples, title="Peptides", normalized=False, save=BASE)
    
    if venn_type == "protein_cluster":
        protein_clusters = defaultdict(set)
        grouped_loci_filt = utils.load(os.path.join(BASE,"grouped_loci_filt.pkl.gz"))
        for cluster in grouped_loci_filt:
            for sample in cluster.quantification.keys():
                protein_clusters[sample].add(cluster.cluster_id)
        venn.make_vens(protein_clusters, samples, title="Protein_Clusters", normalized=False, save=BASE)

def do_pca(args):
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt
    from metaproteomics.analysis import build_loci
    
    metadata = build_loci.read_metadata(os.path.join(BASE, args.metadata))
    df = pd.read_csv(os.path.join(BASE, "X.csv"), index_col=0)
    
    GROUP = args.colorby
    PLOT_PCS = tuple(map(int,args.plot_pcs.split(',')))
    NUM_PCS = max(PLOT_PCS)
    X_scaled = pd.DataFrame(preprocessing.scale(df), index=df.index, columns = df.columns)
    pca = PCA(n_components=NUM_PCS)
    X_r = pd.DataFrame(pca.fit_transform(X_scaled), index=df.index, columns = list(range(1,NUM_PCS+1)))
    print('explained variance ratios: %s' % str(pca.explained_variance_ratio_))
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    groups = np.unique(list(metadata.loc[GROUP]))
    num_colors = len(groups)
    cm = plt.get_cmap("Set1")
    colors = [cm(1.*i/num_colors) for i in range(num_colors)]
    
    for i,group in enumerate(groups):
        samples = list(metadata.columns[metadata.loc[GROUP]==group])
        data = X_r[X_r.index.isin(samples)]
        ax.scatter(data[PLOT_PCS[0]], data[PLOT_PCS[1]], c=[colors[i]]*len(data), label=group, s=50)
    
    
    plt.legend(loc='best',title=GROUP)
    plt.title('PCA of Protein Clusters')
    plt.xlabel("PC{0}: {1:.2f}%".format(PLOT_PCS[0], pca.explained_variance_ratio_[PLOT_PCS[0]-1]*100))
    plt.ylabel("PC{0}: {1:.2f}%".format(PLOT_PCS[1], pca.explained_variance_ratio_[PLOT_PCS[1]-1]*100))
    fig.savefig(os.path.join(BASE, args.out + ".png"))
    fig.savefig(os.path.join(BASE, args.out + ".pdf"))


def main(args):
    metadata = build_loci.read_metadata(os.path.join(BASE, args.metadata))
    
    if args.run_locally:
        run_pipeline(metadata, params, quiet = args.quiet)
    else:
        make_analysis_job(BASE, params['queue'])

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-locally', help="Don't submit garibaldi job. Launch locally", action='store_true')
    parser.add_argument('--quiet', help="Don't display progress bars", action='store_true')
    parser.add_argument('-m','--metadata', help="specify metadata.csv", default="metadata.csv")
    
    """        
    subparsers = parser.add_subparsers(help='Further analysis', dest='analysis')
    subparsers.required = False
    pca_parser = subparsers.add_parser("pca")
    pca_parser.add_argument("--plot-pcs", help="Which PCs to plot. Ex: '1,2'", default="1,2")
    pca_parser.add_argument("--colorby", help="Color samples by group `group`", required=True)
    pca_parser.add_argument('-o',"--out", help="filename to save plot as", default = "pca")   

    venn_parser = subparsers.add_parser("venn")
    venn_parser.add_argument("--samples", help="Samples to use. Ex: `BioGly1,Biogly2'", required=True)
    venn_parser.add_argument("--type", choices = ['peptide','protein_cluster'], help="Data to use in venn", required=True)
    """

    args = parser.parse_args()
    print(args)
    
    """
    # Run the requested pipeline
    switcher = {"pca": do_pca, "venn": do_venn, None: main}
    if args.analysis not in switcher:
        raise ValueError("Incorrect sub-program: '{}'".format(args.analysis))
    switcher[args.analysis](args)
    
    """
    main(args)
