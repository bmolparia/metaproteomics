#!/usr/bin/env python3
"""
This will take all subdirectories within a given directory,
search them for DTASelect-filter.txt files,
parse the DTASelect-filter.txt file for forward loci,
then lookup the GenID for the loci.

For each subdirectory (coresponding to a mass spec run) a column will 
be created in the output file, with rows corresponding to GenIDs
and spectoral counts for each gene will be stored.

The output is a tab-delimited csv that the GSEA program can
read as input.
"""

import os
import pandas as pd
from sys import argv
from pymongo import MongoClient
from metaproteomics.file_processing import blazmass_tools

#setup mongo database handling
mongo_host = 'wl-cmadmin'
mongo_port = 27018
seqdb_name = 'SeqDB_UniprotHuman_HT_2015_02'
seqdb_coll = 'SeqDB_UniprotHuman_HT_2015_02'
protdb_name = 'ProtDB_UniprotHuman_HT_2015_02'
protdb_coll = 'ProtDB_UniprotHuman_HT_2015_02'

seqDB = MongoClient(mongo_host, mongo_port)[seqdb_name][seqdb_coll]
protDB = MongoClient(mongo_host, mongo_port)[protdb_name][protdb_coll]

def main():
    head_path = os.getcwd()  
    tail_path = 'DTASelect-filter.txt'
    out_path = os.path.abspath('out.txt')

    # Argument Handling    
    if len(argv) == 1: #Retain default input and output values
        print('No path entered.  Using current directory')
    if len(argv) > 1:  #First arugment is input value
        head_path = os.path.abspath(os.path.expanduser(argv[1]))
    if len(argv) > 2:  #Second argument is output value
        out_path = os.path.abspath(os.path.expanduser(argv[2]))
    if len(argv) > 3:  #Too many arguments!
        raise ValueError('Too many arguments given')
    
    print('Finding DTASelect-filter.txt files in subdirectories of %s' % head_path)
    print('Output will be saved to %s' % out_path)

    #Find sudirectories and store as list
    runs = [directory for directory in sorted(os.listdir(head_path)) if 
        os.path.isdir(os.path.join(head_path, directory))]


    data_frames = []
    
    print('Reading...')
    # get the data from the files
    gen_names = dict()     
    for i, run in enumerate(runs):
        # set path for current run and initalize parser        
        path = os.path.join(head_path, run, tail_path)        
        parser = blazmass_tools.dta_select_parser(path)
        
        d = dict()
        
        try:        
            n = next(parser)
        except(StopIteration):
            next

        # Loop through each record while there continue to be more
        while (n):
            for locus in n['loci']:  # Look at all loci within a record
                if locus['reverse']:  # Skip resverse Loci
                    next
                record = protDB.find_one({'_id' : locus['Locus']})
                if record['d'].find('GN=') != -1:  
                    #GeneID is stored in Def Line ('d') after GN= and before a space
                    gen_id = record['d'].split('GN=')[1].split(' ')[0]
                    if gen_id in d:
                        d[gen_id] += locus['Spectrum Count']
                    else:
                        d.update({gen_id : locus['Spectrum Count']})
                    
                    # Get descrption if needed                    
                    if gen_id not in gen_names:
                        #descritption stored between | and OS=
                        gen_name = record['d'].split('|')[-1].split('OS=')[0]
                        gen_names.update({gen_id : gen_name})
            
            try: #Final record will throw error when trying to advance beyond
                n = next(parser)
            except(StopIteration):
                break
        #store the dictinary as a pandas dataframe
        data_frames.append(pd.DataFrame.from_dict(d, orient='index'))
        data_frames[i].columns = [run]
        print('Finished reading %s! Scanning next run...' % run)
    
    #put the dataframes together, insert descriptions, and write to file
    result = pd.concat(data_frames, axis = 1)
    result.insert(0,'DESCRIPTION', pd.DataFrame.from_dict(gen_names, orient='index'))
    print('Writing to file \'%s\'...' % out_path)
    open(out_path, 'w').close()
    result.fillna(0).to_csv(out_path, index_label = 'NAME', sep = '\t')

    print("Done!")

if __name__ == '__main__': 
    main()