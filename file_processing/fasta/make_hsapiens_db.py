#!/usr/bin/env python3

# Generates a new FASTA file -- a subset of IndexDB ProtDB matching regex "Homo sapiens.*"
# 	- including "Reverse_" proteins)
#	- maintains ProtDB_072114 protein ID numbers
#
# Sandip Chatterjee
# 11/2/14

from pymongo import MongoClient

def main():
	client = MongoClient('localhost', 27018)
	db = client['ProtDB_072114']
	coll = db['ProtDB_072114']

	count = 0
	# for protDB_record in coll.find({"o":{'$regex':"/Bacteroides fragilis.*/i"}}):
	with open('071414_indexDB_subset_HsapiensDB_reversed.fasta', 'w') as f:
		for protDB_record in coll.find({"o":{'$regex':"Homo sapiens.*"}}):
			count += 1
			f.write('>'+str(protDB_record['_id'])+'||'+protDB_record['d']+'|'+protDB_record['r']+'|[]'+'\n') #ignoring any GO terms associated with ProtDB record
			f.write('\n'.join(split_string_by_n(protDB_record['s'],80))+'\n')

	print(count, "protein records retrieved")

def split_string_by_n(long_string,n):
	while long_string:
		yield long_string[:n]
		long_string = long_string[n:]

if __name__ == '__main__':
	main()