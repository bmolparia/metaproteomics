# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:57:38 2015

@author: gstupp
"""

import os
from collections import defaultdict
#folder_path = "/home/gstupp/projects/Wolan/ana_probe_cmkGLP/2015_04_15_compil_subset/kC_1mil_3"

def parse_kClust(folder_path):
    header_path = os.path.join(folder_path, 'headers.dmp')
    clusters_path = os.path.join(folder_path, 'clusters.dmp')
    
    # read kClust# to protID mapping
    with open(header_path) as file:
        headers = dict()
        for line in file:
            prot_id = int(line.split('>',1)[1].split('||',1)[0])
            kClust_id = int(line.split('>',1)[0].strip())
            headers[kClust_id] = prot_id
    # Read clusters file
    with open(clusters_path) as file:
        groups = defaultdict(list)
        next(file) #header
        for line in file:
            prot_id_k, group_id_k = [int(x) for x in line.strip().split()]
            prot_id = headers[prot_id_k]
            group_id = headers[group_id_k]
            groups[group_id].append(prot_id)
    return dict(groups)

def parse_kClust_compil(folder_path):
    # format as list of dicts, like so:
    # {'_id': 123, 'pID': [123,133,353]}
    groups = parse_kClust(folder_path)
    json_groups = []
    for _id, group in groups.items():
        json_groups.append({'_id': _id, 'pID': group})
    return json_groups


'''
from itertools import chain
inv_groups = {v:k for (k,v) in chain(*[list(zip([k]*len(v), v)) for (k,v) in groups.items()])}
'''