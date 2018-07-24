# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:20:03 2015

@author: gstupp

Pull out parent ion intensities for each scan from ms2 files.
Pull out scan intensities for each peptide from DTASelect file
Parse a DTASelect-filter.txt file. Grab top N scans for each peptide.
Fill in ion quantification for each peptide.

According to this paper: http://www.nature.com/nbt/journal/v28/n1/full/nbt.1592.html
this is something like the AUC. Although the AUC is really quantified from the elution 
chromatogram. Ideally we, would want the intensities of the fragment ions, but
those are nowhere to be found...

Input:
    ms2_path = path to folder containing ms2 files. Globs *.ms2
    dtaselect_path = path to folder containing DTASelect.txt and DTASelect-filter.txt
    out_file_path = path to output DTASelect-filter.txt (default overwrite orig)
"""
import os
import glob
from itertools import chain
from metaproteomics.utils import get_lcstep
from collections import defaultdict
from file_processing import blazmass_tools

#ms2_path = '/mongoc/gstupp/DTASelect/01_2015_mass_spec/H1_11082014/1108_Phe2_2014_12_12_10_29176/'
#dtaselect_path = '/mongoc/gstupp/DTASelect/01_2015_mass_spec/H1_11082014/1108_Phe2_2014_12_12_10_29176/indexDB_search_10ppm_50ppmfrag/dtaselect_results_sfp0.01_p2'

def parse_ms2(ms2_path):
    # store dict of step.scan.intensity
    # assumes ms2 file name is something like "12112014_H1_1108_Gly1_2.ms2"
    # where the step is everything between the last '_' and the last '.'
    ms2_paths = glob.glob(os.path.join(ms2_path,'*.ms2'))
    scan_int = dict()
    for ms2_path in ms2_paths:
        step = get_lcstep(os.path.basename(ms2_path))
        scan_int[step] = dict()
        
        ms2 = open(ms2_path)
        ms2 = filter(lambda line: line.startswith('S') or line.count('PrecursorInt'), ms2)
        for line in ms2:
            scan = int(line.strip().split()[1])
            line = next(ms2)
            intensity = float(line.strip().split()[2])
            scan_int[step][scan] = intensity
    return scan_int

def parse_dtaselect(dtaselectuniq_path, scan_int):
    '''
    dtaselectuniq_path:
    12102014_H1_1108_Phe2_7.159.159.2	0.23381054	K.LKKGVLIAAIIRK.G
    12102014_H1_1108_Phe2_9.159.159.2	0.2337929	3   K.LKKGVLIAAIIRK.G
    ...
    
    grep $'^D\t' DTASelect.txt | cut -f 2,6,13 | sort | uniq > DTASelect.txt.uniq
    '''
    dta = open(dtaselectuniq_path)
    pep_step_scan = defaultdict(list)
    for line in dta:
        line = line.strip().split()
        peptide = line[2].split('.')[1]
        name = line[0]
        step = get_lcstep(name)
        scan = int(name.split('_')[-1].split('.')[1])
        if peptide in pep_step_scan:
            if any(x[2] == step and x[3] == scan for x in pep_step_scan[peptide]):
                continue
        try:
            intensity = scan_int[step][scan]
        except KeyError:
            print('Failed to find: ' + str(step) + ' ' + str(scan) + ' ' + peptide)
            continue
        score = float(line[1])
        pep_step_scan[peptide].append((score, intensity, step, scan))
    pep_step_scan = dict(pep_step_scan)
    for k,pep_entry in pep_step_scan.items():
        pep_entry.sort(key = lambda x: x[0], reverse = True)
        pep_step_scan[k] = tuple(x[1] for x in pep_entry)
    return pep_step_scan

#sum(pep_step_scan['ALLDQLSQIVADQLENGGEITLPGVGK'][:9])

def dta_select_ion_quant(dtaselectfilter_path, pep_step_scan):
    '''
    Parse a DTASelect-filter.txt file. Fill in ion quantification for each peptide.
    '''    
    out = []
    with open(dtaselectfilter_path) as f:
        for line in f:
            if (line.startswith('\t') or line.startswith('*')) and line.count('Proteins') == 0:
                line = line.split()
                pep_seq = line[-1].split('.')[1]
                redundancy = int(line[-2])
                intensity = sum(pep_step_scan[pep_seq][:redundancy])
                if line[0] == '*':
                    line[7] = '%.1f' % intensity
                    out.append('\t'.join(line) + '\n')
                else:
                    line[6] = str(intensity)
                    out.append('\t' + '\t'.join(line) + '\n')
            else:
                out.append(line)
    return out

def main(ms2_path, dtaselect_path, out_file_path = None):
    dtaselectfilter_path = os.path.join(dtaselect_path, 'DTASelect-filter.txt')
    dtaselecttxt_path = os.path.join(dtaselect_path, 'DTASelect.txt')
    dtaselectuniq_path = dtaselecttxt_path + ".uniq"

    # Parse ms2 files
    scan_int = parse_ms2(ms2_path)
        
    # Make a much smaller DTASelect.txt containing:
    # the scans each peptide was ID'd in along with their score
    os.system("grep -P '^D\t' " + dtaselecttxt_path + " | cut -f 2,6,13 | sort -k 1 | uniq > " + dtaselectuniq_path)
    pep_step_scan = parse_dtaselect(dtaselectuniq_path, scan_int)
    
    new_dtaselect_filter = dta_select_ion_quant(dtaselectfilter_path, pep_step_scan)
    
    if out_file_path == None:
        out_file_path = os.path.join(dtaselect_path, 'DTASelect-filter.txt')
    if out_file_path.count(os.sep) == 0:
        out_file_path = os.path.join(dtaselect_path, out_file_path)
    out_file_path = os.path.abspath(out_file_path)
    with open(out_file_path,'w') as out_file:
        out_file.writelines(new_dtaselect_filter)
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('ms2_path', type = str, help = 'path to folder containing ms2 files. Globs *.ms2')
    parser.add_argument('dtaselect_path', type = str, help = 'path to folder containing DTASelect.txt and DTASelect-filter.txt')
    parser.add_argument('-o','--out_file_path', type = str, help = 'path to output DTASelect-filter.txt (default overwrite orig)')
    args = parser.parse_args()
    main(os.path.abspath(args.ms2_path), os.path.abspath(args.dtaselect_path), args.out_file_path)
    
    
    
    