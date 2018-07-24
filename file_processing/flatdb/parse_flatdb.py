#!/usr/bin/env python

# parse_flatdb.py
# v1, 6/5/14
# Sandip Chatterjee

import sys
import multiprocessing as mp
import json
from time import time

def main():
	flatdb_file = sys.argv[1]
	split_key = sys.argv[2]
	chunk_size = int(sys.argv[3])

	if split_key == 'mass':
		split_position = 0
	elif split_key == 'sequence':
		split_position = 1
	else:
		print "Invalid flatdb split key -- exiting"
		sys.exit()

	pool = mp.Pool()
	
	start_time = time()
	if split_key == 'mass':
		JSON_objects = pool.imap(process_chunk_massdb, read_flatfile(flatdb_file,split_position),chunk_size)
	elif split_key == 'sequence':
		JSON_objects = pool.imap(process_chunk_seqdb, read_flatfile(flatdb_file,split_position),chunk_size)

	print "Time taken to make imap generator:", time()-start_time, "seconds"

	start_time = time()
	print "Saving to new JSON file", flatdb_file.replace('.flatdb','.json')
	with open(flatdb_file.replace('.flatdb','.json'),'wb') as f:
		for obj in JSON_objects:
			f.write(json.dumps(obj)+'\n')
		print "About to close multiprocessing pool. total time taken:", time()-start_time, "seconds"
		pool.close()
		pool.join()

	print "Finished"



def read_flatfile(DB_file,split_position):

	current_chunk = []

	with open(DB_file,'rb') as f:

		line = f.readline()
		while True:
			current_chunk = []
			if not line:
				break
			else:
				current_peptide = line.split()[split_position]
			while True:
				current_chunk.append(line)
				line = f.readline()
				if not line or line.split()[split_position] != current_peptide:
					break

			yield current_chunk

def process_chunk_massdb(current_chunk):

	new_JSON_obj = {}
	s = []
	new_JSON_obj['_id'] = int(current_chunk[0].split()[0])
	for line in current_chunk:
		s.append(line.split()[1])
	new_JSON_obj['s'] = list(set(s))
	return new_JSON_obj

def process_chunk_seqdb(current_chunk):

	new_JSON_obj = {}
	p = []
	new_JSON_obj['_id'] = current_chunk[0].split()[1]

	for line in current_chunk:
		new_JSON_parent_obj = {}
		line_split = line.split()
		new_JSON_parent_obj['l'] = line_split[2]
		new_JSON_parent_obj['r'] = line_split[3]
		new_JSON_parent_obj['i'] = int(line_split[4].lstrip('r'))
		new_JSON_parent_obj['o'] = int(line_split[5])
		if 'r' in line_split[4]:
			new_JSON_parent_obj['d'] = True
		p.append(new_JSON_parent_obj)
	new_JSON_obj['p'] = p
	return new_JSON_obj

if __name__ == '__main__':
	main()
