"""
Some boilerplate venn diagram functions
"""
from matplotlib_venn import venn3, venn2, _venn3, _venn2
import pylab
import numpy as np
import os

def venn(*args, **kwargs):
    # cheap wrapper for matplotlib_venn.venn2/venn3 that calls the appropriate one
    x = args[0] if len(args) == 1 else args[1]
    if len(x) == 2:
        venn2(*args, **kwargs)
    elif len(x) == 3:
        venn3(*args, **kwargs)
    else:
        raise ValueError("Only supports 2 or 3 items")

def make_vens(sample_dict, samples, title='', save=None, normalized=False):
    # sample_dict: {sample: values_to_use}
    # samples: list of sample to use for venn. Is a subset of sample_dict.keys()
    # title: title for plot
    # save: folder to save plot in, or None
    # normalized: boolean normalize total to 100

    if len(samples) not in [2,3]:
        print("2 or 3 samples required")
        return None
    pylab.figure(figsize=(6,6))
    ax = pylab.axes()
    p,l = zip(*[(dta,sample) for sample,dta in sample_dict.items() if sample in samples])
    n_str = ""
    if normalized:
        n_str = "_norm"
        p = np.array(_venn3.compute_venn3_subsets(*p)) if len(samples) == 3 else np.array(_venn2.compute_venn2_subsets(*p))
        p = (100*p/sum(p)).round()
    venn(p,l, ax=ax)
    pylab.title(title)
    if save:
        file_name = ','.join(samples) + "_" + title + n_str + ".pdf"
        pylab.savefig(os.path.join(save,file_name), bbox_inches='tight')

def boolean_make_vens(comparisons, sample_dict, join_type, 
                      title='', names = None, save=None):
    #join_type = 'union' or 'intersection'
    f = getattr(set, join_type)
    members = [f(*[v for k,v in sample_dict.items() if k in group]) for group in comparisons]
    pylab.figure(figsize=(6,6))
    ax = pylab.axes()
    if not names:
        names = [join_type + '(' + ','.join(group) + ')' if len(group)>1 else group[0] for group in comparisons]
    venn(members, names, ax=ax)
    pylab.title(title)
    if save:           
        file_name = ','.join(names) + "_" + title + ".pdf"
        pylab.savefig(os.path.join(save,file_name), bbox_inches='tight')

def run_example():
    sample_dict = {'a':{1,2,3,4,5,6}, 'b':{3,4,5,6,7}, 'c':{1,5,6,7,8,9,10,11,12}, 'd':{0,1,2}}
    make_vens(sample_dict, ['a','b'], 'title')
    make_vens(sample_dict, ['a','b'], 'title', save=False, normalized=True)
    boolean_make_vens([['a','b'], ['c','d']], sample_dict, "union", 'title') # union of a and b vs. union of c and d