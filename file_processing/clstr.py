"""
Functions for parsing cd-hit output
parse_clstr: parse clstr file
parse_clstr_compil: parse clstr file for inserting into mongodb. specific to compil format (having '||' as separator and starting with protID)
mongo_clstr: insert clstr file into a mongodb. IPA = True -> calls annotate_group_ipa on each doc
annotate_group_ipa: annotate each cluster with the union of all ipa terms for all protIDs in that group
add_annotations: for adding annotations to an already existing db. iterates through all docs and updates them with the annotations
"""

#clstr_file = "/mongoc/gstupp/cdhit/seen_complil/seen_compil_g1.clstr"
#clstr_file = "/mongoc/clustering/071414_ComPIL_forwardonly_cluster10/071414_ComPIL_forwardonly_cluster10_0.7.clstr"

def cluster_chunk(clstr_file):
    lines = []
    for line in open(clstr_file):
        if line[0] == '>':
            if lines:
                yield lines
                lines = []
        else:
            lines.append(line)
            
def parse_clstr(clstr_file, name = 0, split_at = '||'):
    # split_at use '...' for non compil
    # name == 0: group name is an integer that increments, 1: parent's id (annotated with *)
    c = cluster_chunk(clstr_file)
    groups = dict()
    for idx,chunk in enumerate(c):
        if name == 0:
            root = idx
        elif name == 1:
            root = [line.split('>')[1].split(split_at)[0] for line in chunk if line.rstrip().endswith('*')][0]
        else:
            print('bad option')
            return
        ids = [line.split('>')[1].split(split_at)[0] for line in chunk]   
        if split_at == '||':
            ids = [int(id_) for id_ in ids]
            root = int(root)
        groups[root] = ids
    return groups

def parse_clstr_compil(clstr_file):
    c = cluster_chunk(clstr_file)
    for chunk in c:
        root = int([line.split('>')[1].split('||')[0] for line in chunk if line.rstrip().endswith('*')][0])
        ids = [int(line.split('>')[1].split('||')[0]) for line in chunk]
        yield {'_id': root, 'pID': ids}


def mongo_clstr(clstr_file, name, host = 'wl-cmadmin', port = None, ipa = True):
    #name = "071414_ComPIL_forwardonly_0_7"
    from pymongo import MongoClient
    coll = MongoClient(host = host, port = port)[name][name]
    
    if ipa:
        coll.insert(map(annotate_group_ipa, parse_clstr_compil(clstr_file)))
    else:
        coll.insert(parse_clstr_compil(clstr_file))
        
    coll.ensure_index('pID')

def annotate_group_ipa(group):
    from analysis import functional_analysis
    functional_analysis.init()
    ipa = functional_analysis.get_annotations_from_protIDs(group['pID'], return_ipa = True)
    group['ipa'] = list(ipa) if ipa else []
    return group

def add_annotations(name, host = 'wl-cmadmin', port = None):
    # for adding annotations to an already existing db
    from pymongo import MongoClient
    coll = MongoClient(host = host, port = port)[name][name]
    
    for doc in coll.find({'ipa': {'$exists': False}}):
        doc = annotate_group_ipa(doc)
        coll.update({'_id': doc['_id']}, doc, w = 0)

#pa = map(annotate_group_ipa, p)

#%% 
if 0:
    '''
    #%%
    #from itertools import chain
    #groups = parse_clstr(clstr_file, name = 1)
    #inv_groups = {v:k for (k,v) in chain(*[list(zip([k]*len(v), v)) for (k,v) in groups.items()])}
    #%%
    g1 = [inv_groups[p] for p in ps[4]['root_protID']]
    g2 = [inv_groups[p] for p in ps[6]['root_protID']]
    m1 = [groups[x] for x in set(g1)]
    m2 = [groups[x] for x in set(g2)]
    #%%
    from analysis import build_loci_from_all_peptides
    for idx1,p1 in enumerate(ps):
        for idx2,p2 in enumerate(ps):
            if p1 != p2:
                s = build_loci_from_all_peptides.similarity(p1['peptide_seq'], p2['peptide_seq'])
                if s > .5:
                    print(str(idx1) + ',' + str(idx2))
                    break
    4,6
    6,4
    7,20
    8,17
    9,4
    10,24
    11,15
    12,4
    14,18
    15,11
    17,8
    18,14
    19,23
    20,7
    '''
