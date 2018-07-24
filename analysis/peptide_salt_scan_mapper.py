"""
Created on Fri Mar 20 13:20:03 2015

@author: gstupp

Pull out score_salt-step_scan for each peptide from DTASelect file
Parse a DTASelect-filter.txt file. Grab top N scans for each peptide.

Input:
    dtaselect_path = path to folder containing DTASelect.txt and DTASelect-filter.txt
"""
import os
from collections import defaultdict
from itertools import chain
from collections import Counter

#dtaselect_path = '/mongoc/gstupp/DTASelect/01_2015_mass_spec/H1_11082014/1108_Phe2_2014_12_12_10_29176/indexDB_search_10ppm_50ppmfrag/dtaselect_results_sfp0.01_p2'
try:
    get_ipython().magic('pylab inline')
except:
    pass
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as fm
fpath = '/usr/share/fonts/Roboto/Roboto-Regular.ttf'
#fpath = '/mongoa/ipynb/sandip/HelveticaNeue.ttf'
prop36 = fm.FontProperties(fname=fpath, size=36)
prop32 = fm.FontProperties(fname=fpath, size=32)
prop30 = fm.FontProperties(fname=fpath, size=30)
prop24 = fm.FontProperties(fname=fpath, size=24)
prop20 = fm.FontProperties(fname=fpath, size=20)
prop18 = fm.FontProperties(fname=fpath, size=18)
prop16 = fm.FontProperties(fname=fpath, size=16)
prop14 = fm.FontProperties(fname=fpath, size=14)

color = '#3498DB'
opacity = 1
file_extension = '.pdf' # .pdf, .jpg, .png


def parse_dtaselect(dtaselectuniq_path):

    dta = open(dtaselectuniq_path)
    pep_step_scan = defaultdict(list)
    for line in dta:
        line = line.strip().split()
        peptide = line[2].split('.')[1]
        name = line[0]
        step = int(name.split('_')[-1].split('.')[0])
        scan = int(name.split('_')[-1].split('.')[1])
        charge = int(name.split('_')[-1].split('.')[3])
        score = float(line[1])
        pep_step_scan[peptide].append((score, step, scan, charge))
    pep_step_scan = dict(pep_step_scan)
    for k,pep_entry in pep_step_scan.items():
        pep_entry.sort(key = lambda x: x[0], reverse = True)
    return pep_step_scan
    
    '''
    Example output:
    pep_step_scan['ADPNIPDVQVIVK']
    Out[96]: 
    [(0.39606094, 5, 6094, 2),
     (0.38992006, 6, 5664, 2),
     (0.38867128, 5, 6070, 2),
     (0.36590993, 6, 5714, 2),
     (0.34838128, 6, 5871, 2)]
     '''

def dta_select_ion_quant(dtaselectfilter_path, pep_step_scan):
    '''
    Parse a DTASelect-filter.txt file. 
    Filter pep_step_scan to keep only the top X, specified by its redundancy
    '''    
    pep_step_scan_filtered = dict()
    with open(dtaselectfilter_path) as f:
        for line in f:
            if (line.startswith('\t') or line.startswith('*')) and line.count('Proteins') == 0:
                line = line.split()
                pep_seq = line[-1].split('.')[1]
                redundancy = int(line[-2])
                pep_step_scan_filtered[pep_seq] = pep_step_scan[pep_seq][:redundancy]
                
    return pep_step_scan_filtered

def main(dtaselect_path):
    dtaselectfilter_path = os.path.join(dtaselect_path, 'DTASelect-filter.txt')
    dtaselecttxt_path = os.path.join(dtaselect_path, 'DTASelect.txt')
    dtaselectuniq_path = dtaselecttxt_path + ".uniq"

    # Make a much smaller DTASelect.txt containing:
    # the scans each peptide was ID'd in along with their score
    os.system("grep -P '^D\t' " + dtaselecttxt_path + " | cut -f 2,6,13 | sort -k 1 | uniq > " + dtaselectuniq_path)
    pep_step_scan = parse_dtaselect(dtaselectuniq_path)
    
    pep_step_scan_filtered = dta_select_ion_quant(dtaselectfilter_path, pep_step_scan)
    
    psm_ids_set = get_distinct_psm_ids(pep_step_scan_filtered)
    
    histogram = make_LCStep_histogram(psm_ids_set)
    
    return histogram

def get_distinct_psm_ids(pep_step_scan_filtered):
    return {'_'.join(map(str,x[1:])) for x in chain(*pep_step_scan_filtered.values())}
        
def make_LCStep_histogram(psm_ids_set):
    full_lcstep_count = []
    for psm in psm_ids_set:
        full_lcstep_count.append(psm.split('_')[0])

    return Counter(full_lcstep_count)

def plot_salt(histogram, name, ymaxlimit=None):

    fig, ax = plt.subplots()
    fig.set_size_inches(6,5)
#     labels, values = zip(*sorted(hist.items()))
    labels, values = zip(*sorted(histogram.items(), key=lambda x: int(x[0])))

    indexes = np.arange(len(labels))
    width = 1

    plt.bar(indexes, values, width, color=color, alpha=opacity)

    fmt = plt.ScalarFormatter(useOffset=False)
    plt.gca().xaxis.set_major_formatter(fmt)
    ax.set_xlabel('LC Step #', fontproperties=prop32)
    ax.set_ylabel('# high quality PSMs', fontproperties=prop32)
    if ymaxlimit:
        ax.set_ylim(0,ymaxlimit)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop24)

    plt.title('Peptides identified per Liquid Chromatography Step\n{}'.format(name),fontproperties=prop24, y=1.1)
    ax.set_xticks(indexes+width/2)
    ax.set_xticklabels(labels,fontproperties=prop24)
    #plt.tight_layout()
    plt.show()
