# -*- coding: utf-8 -*-
"""
Functions for calculating the SpI of a protein from spectral counts in two groups
and calculating a p-value for the change by permutation tests

"""
import numpy as np
def perm_test_SpI(group1, group2, nRep = 10000):
    """Permutation test to calc p-value for difference between group1 and group2 using SpI"""
    observedSpI = SpI(group1, group2)
    combined = np.concatenate([group1, group2])
    num1 = len(group1)
    resampleSpI = np.zeros(nRep, dtype='float')
    for idx in range(nRep):
        permutedCombined = np.random.permutation(combined)
        resampleSpI[idx] = SpI(permutedCombined[:num1], permutedCombined[num1:])
        
    pVal = (np.sum(resampleSpI > observedSpI) + np.sum(resampleSpI < -observedSpI))/float(nRep)
    return observedSpI, resampleSpI, pVal

def perm_test_median(group1, group2, nRep = 10000):
    """Permutation test to calc p-value for difference between group1 and group2 using difference between medians"""
    group1_0 = [0 if x == None else x for x in group1]
    group2_0 = [0 if x == None else x for x in group2]
    observed_median_diff = np.median(group1_0) - np.median(group2_0)
    combined = np.concatenate([group1_0, group2_0])
    num1 = len(group1)
    resample_md = np.zeros(nRep, dtype='float')
    for idx in range(nRep):
        permutedCombined = np.random.permutation(combined)
        resample_md[idx] = np.median(permutedCombined[:num1]) - np.median(permutedCombined[num1:])
    
    return observed_median_diff, resample_md

def perm_test_mean(group1, group2, nRep = 10000):
    """Permutation test to calc p-value for difference between group1 and group2 using difference between means"""
    group1_0 = [0 if x == None else x for x in group1]
    group2_0 = [0 if x == None else x for x in group2]
    observed_median_diff = np.nanmean(group1_0) - np.nanmean(group2_0)
    combined = np.concatenate([group1_0, group2_0])
    num1 = len(group1)
    resample_md = np.zeros(nRep, dtype='float')
    for idx in range(nRep):
        permutedCombined = np.random.permutation(combined)
        resample_md[idx] = np.nanmean(permutedCombined[:num1]) - np.nanmean(permutedCombined[num1:])
    
    return observed_median_diff, resample_md

def SpI(group1, group2):
    # http://pubs.acs.org/doi/full/10.1021/pr070271%2B
    # group1 = [3,12,8,None,None,None]
    # group2 = [30,10]
    group1NA = [np.NaN if x == None or x == 0 else x for x in group1]
    group2NA = [np.NaN if x == None or x == 0 else x for x in group2]
    num_det1 = len([x for x in group1 if x])
    num_det2 = len([x for x in group2 if x])
    if num_det1 == 0:
        return -1.0
    elif num_det2 == 0:
        return 1.0
    msc = np.nansum([np.nanmean(group1NA), np.nanmean(group2NA)])
    left = ( np.nanmean(group1NA) / msc ) * (num_det1/len(group1))
    right = ( np.nanmean(group2NA) / msc ) * (num_det2/len(group2))
    return left - right