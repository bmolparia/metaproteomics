





# Normalization refs
# http://proteomesci.biomedcentral.com/articles/10.1186/1477-5956-11-S1-S13
# http://www.ncbi.nlm.nih.gov/pubmed/23808607
# http://online.liebertpub.com/doi/full/10.1089/omi.2013.0010


if False:
    #%% boxplot of peptide quant
    import numpy as np
    import matplotlib.pyplot as plt
    sp = {sample.sample_name:np.array(list(sample.pep_quant.values())) for sample in samples.values()}
    sp_log = {sample:np.log(v+1) for sample,v in sp.items()}
    plt.boxplot(list(sp_log.values()), labels=list(sp_log.keys()))
    {sample.sample_name:np.sum(np.array(list(sample.pep_quant.values()))) for sample in samples.values()}
    #%% boxplot of grouped_loci quant
    gl = dict()
    for locus in grouped_loci:
        for key,value in locus.quantification.items():
            gl.setdefault(key,[]).append(value)
    gl_log = {sample:np.log(np.array(v)+1) for sample,v in gl.items()}
    plt.boxplot(list(gl_log.values()), labels=list(gl_log.keys()))
    
    {sample:sum(v) for sample,v in gl.items()}
    

#%% operations on matrix formed from grouped_loci
if False:
    
    # distance matrix
    import scipy.spatial.distance as ssd
    D = ssd.squareform(ssd.pdist(X, metric = "correlation"))
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    im = ax.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    ax.set_xticks(list(range(18)))
    ax.set_xticklabels(list(df.index), rotation=90)
    ax.set_yticks(list(range(18)))
    ax.set_yticklabels(list(df.index), rotation=0)

#%%





