"""
Created on Tue Jan 20 10:50:34 2015

@author: gstupp

Figuring out how to annotate loci
Work in progress


From home run:
$ ssh -fNL 27028:localhost:27018 gbwl
$ ssh -fNL 27027:localhost:27017 gbwl
where gbwl is username@garibaldi.scripps.edu -W wl-cmadmin.scripps.edu:22

"""

#HOST = 'localhost' 
HOST = 'wl-cmadmin'

from pymongo import MongoClient
from itertools import chain
from collections import Counter
import taxonomy
t = taxonomy.Taxonomy(host = HOST, port=27017)

client = MongoClient(HOST, 27017)
taxDB = client.TaxDB_072114.TaxDB_072114
hashDB = client.HashDB_072114.HashDB_072114
domainDB = client.DomainDB_072114.DomainDB_072114

protDB = MongoClient(HOST, 27018).ProtDB_072114.ProtDB_072114

from file_processing import blazmass_tools
from file_processing.dta import build_proteins_from_peptides

def get_domains(p):
    hashes = [x['_id'] for x in hashDB.find({'pID':{'$in': p['forward_loci']}}, {'_id': True })]
    domain_result = domainDB.find({'_id':{'$in': hashes}})
    domain_dict = dict([(x['_id'],x['d']) for x in domain_result]) # dictionary where hashes are keys
    domain_list = list(chain(*[x for x in chain(*[x.values() for x in domain_dict.values()])])) # flat list of all domains in all hashes

    #### Now same thing just for parent loci
    parent_hashes = [x['_id'] for x in hashDB.find({'pID':{'$in': p['parent_forward_loci']}}, {'_id': True })]
    parent_domain_result = domainDB.find({'_id':{'$in': parent_hashes}})
    parent_domain_dict = dict([(x['_id'],x['d']) for x in parent_domain_result])
    parent_domain_list = list(chain(*[x for x in chain(*[x.values() for x in parent_domain_dict.values()])]))

    # get set of Go terms for each hash
    hash_go_dict = {hash:set(chain(*[x['g'] for x in chain(*[x for x in values.values()]) if 'g' in x])) for hash, values in domain_dict.items()}
    go_set = set.intersection(*list(hash_go_dict.values())) if hash_go_dict else {}
    hash_go_dict_counter = Counter([frozenset(x) for x in hash_go_dict.values()]) # What are the unique sets of go terms?
    parent_hash_go_dict = {hash:set(chain(*[x['g'] for x in chain(*[x for x in values.values()]) if 'g' in x])) for hash, values in parent_domain_dict.items()}
    parent_go_set = set.intersection(*list(parent_hash_go_dict.values())) if parent_hash_go_dict else {}
    hash_go_dict_counter = Counter([frozenset(x) for x in parent_hash_go_dict.values()])
    
    p['go_set'] = go_set
    p['parent_go_set'] = parent_go_set
    p['domain_list'] = domain_list
    p['parent_domain_list'] = parent_domain_list
    p['domain_result'] = domain_result

# Get children of protease go term 
def get_go_protease():
    from goatools.obo_parser import GODag
    g = GODag('/home/gstupp/go/go-basic.obo')
    go_term = g.query_term('GO:0008233')
    go_protease = go_term.get_all_children()
    return go_protease

#go_protease = get_go_protease()

dta_select = '/home/gstupp/01_2015_mass_spec/H1_11082014/1108_Gly1_2014_12_15_15_29205/dtaselect_results_sfp0.01_p2/DTASelect-filter.txt'
dta_select = '/home/gstupp/01_2015_mass_spec/102214_HEK293_HCD35/indexDB_search_noProtDB/DTASelect-filter.txt'
ps = build_proteins_from_peptides.main(dta_select, mongo_host=HOST, mongo_port=27018)

ps = [p for p in ps if not p['reverse']]
for idx,p in enumerate(ps):
    p = get_domains(p)

# Number of loci in which the go terms are the same between parent and subsets
len([p for p in ps if p['go_set'] == p['parent_go_set'] and p['go_set'] != {}])
# Number in which there are no GOs for either
len([p for p in ps if p['parent_go_set'] == {} and p['go_set'] == {}])
# Number in which Only in subset
len([p for p in ps if p['go_set'] != p['parent_go_set'] and p['parent_go_set'] == {}])
# Only in parent
len([p for p in ps if p['go_set'] != p['parent_go_set'] and p['go_set'] == {}])





####################################3



ps_loci = set(chain(*[p['forward_loci'] for p in ps]))
ps_ploci = set(chain(*[p['parent_forward_loci'] for p in ps]))

loci_getf = set(chain(*[p['forward_loci'] for p in ps_getf]))
loci = set(chain(*[p['forward_loci'] for p in ps]))

for p in ps:
    if p['set_go'] and len(set.intersection(p['set_go'],go_protease)) > 0:
        p['protease'] = True
    else:
        p['protease'] = False

print(sum([p['protease'] for p in ps]))
set_go = list(set(chain(*[p['set_go'] for p in ps if p['set_go']])))
all_peptides = list(set(chain(*[x['peptide_seq'] for x in ps]))) #set of all peptides in the whole file
all_peptides_lr = list(set(chain(*[[x['Sequence'] for x in p['peptides']] for p in ps])))
all_peptides_lr = sorted(all_peptides_lr, key=lambda x:x.split('.')[1])
from itertools import groupby
diff_peptides = []
grouper = groupby(all_peptides_lr, key = lambda x:x.split('.')[1])
for peptide, group in grouper:
    peptides = list(group)
    if len(peptides)>1:
        diff_peptides.append(peptides)
        

################## 'FAAYIQQSNMESNGK'

# Map GoSlim
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim
from collections import Counter
go_dag = GODag('/home/gstupp/goatools/go.obo')
goslim_dag = GODag('/home/gstupp/goatools3/goslim_generic.obo')
goslim_meta = GODag('/home/gstupp/goatools3/goslim_metagenomics.obo')

ana_DTA = '/home/gstupp/01_2015_mass_spec/H1_11082014/1108_Gly1_2014_12_15_15_29205/dtaselect_results_sfp0.01_p2/DTASelect-filter.txt'
parser = blazmass_tools.dta_select_parser(ana_DTA, small = True)
ps = [get_domains(p) for p in parser]
set_go = set(chain(*[p['set_go'] for p in ps if p['set_go'] is not None]))
for p in ps:
    if p['set_go']:
        p['go_slim'] = set(chain(*[mapslim(go_term, go_dag, goslim_meta)[0] for go_term in p['set_go'] if go_term in go_dag]))
    else:
        p['go_slim'] = None
go_slim = Counter(chain(*[p['go_slim'] for p in ps if p['go_slim']]))
labels = {go_term:go_dag.query_term(go_term).name for go_term in go_slim.keys()}
[labels[go] for (go,x) in go_slim.most_common(n=10)]
    
    
    
import plot_tools



cmap = plt.cm.jet
colors = cmap(np.linspace(0., 1., len(go_slim.keys()))) 
explode = [1 if x < 20 else 0 for x in go_slim.values()]

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(111)
patches = ax.pie(list(go_slim.values()), autopct='%1.1f', explode = explode, colors = colors)
ax.set_title('LCA for each locus, Bfrag BfragDB large ppm')
# Put a legend to the left of the current axis
ax.legend(patches[0], labels, loc='center left', bbox_to_anchor=(-.6, .5))

