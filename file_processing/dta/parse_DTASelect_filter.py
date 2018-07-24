#!/usr/bin/env python3

# DTASelect parser v0.2
#
# Sandip Chatterjee
# 12/9/14

import argparse

def main():
	
	parser = argparse.ArgumentParser()
	parser.add_argument('filepath', help='Path to DTASelect-filter.txt file', type=str)
	parser.add_argument('--saveloci', help='Save ProtDB IDs (forward only) to a file, one ProtDB ID per line (DTASelect_forwardloci)', action='store_true')
	parser.add_argument('--outputprefix', help='String to prepend to filename(s) for output files', type=str)
	parser.add_argument('--savepeptides', help='Save peptide matches (forward only) to a file, one peptide sequence per line (DTASelect_forwardpeptides)', action='store_true')

	args = parser.parse_args()

	filepath = args.filepath

	print("Reading DTASelect output file: {}",format(filepath), end='\n\n')

	# file_lines = read_DTASelect_filter_file(filepath)
	header, data, footer = read_DTASelect_filter_file(filepath)

	print("HEADER", header)
	print('-----------------')
	print("FOOTER", footer)
	print('-----------------')
	print("DATA START", data[:2])
	print('-----------------')
	print("DATA END", data[-2:])
	print('--------------')

	parsed_data_chunks = list(retrieve_loci_noProtDB(data)) 	## exhausting generator because file is small
	chunks_with_peptide_matches = [chunk for chunk in parsed_data_chunks if len(chunk) > 1]
	forward_chunks_with_peptide_matches = filter_chunks_to_forward_only(chunks_with_peptide_matches)
	forward_peptide_matches = get_peptides(forward_chunks_with_peptide_matches)

	protDB_IDs_for_chunks_with_peptide_matches = [chunk[0].split('\t')[0] for chunk in chunks_with_peptide_matches]
	forward_protDB_IDs_for_chunks_with_peptide_matches = [int(ID) for ID in protDB_IDs_for_chunks_with_peptide_matches if not ID.startswith('Reverse_')]
	reverse_protDB_IDs_for_chunks_with_peptide_matches = [int(ID.replace('Reverse_','')) for ID in protDB_IDs_for_chunks_with_peptide_matches if ID.startswith('Reverse_')]

	print('Loci with 1+ peptides matches:', len(chunks_with_peptide_matches))
	print('ProtDB IDs for Loci with 1+ peptides matches:', len(protDB_IDs_for_chunks_with_peptide_matches))
	print('Forward Loci with 1+ peptides matches:', len(forward_protDB_IDs_for_chunks_with_peptide_matches))
	print('Reverse Loci with 1+ peptides matches:', len(reverse_protDB_IDs_for_chunks_with_peptide_matches))


	chunks_with_unique_peptide_matches = [chunk for chunk in chunks_with_peptide_matches if [peptide_match for peptide_match in chunk[1:] if peptide_match[0] == '*']]
	print('Loci with 1+ unique peptide matches:', len(chunks_with_unique_peptide_matches))
	print('--------------')

	if args.saveloci:
		if args.outputprefix:
			prefix = args.outputprefix
		else:
			prefix = ''
		with open(prefix+'_'+'DTASelect_forwardloci','w') as f:
			for locus in forward_protDB_IDs_for_chunks_with_peptide_matches:
				f.write(str(locus)+'\n')

		print('ProtDB IDs saved to file {}_DTASelect_forwardloci'.format(prefix))

	if args.savepeptides:
		if args.outputprefix:
			prefix = args.outputprefix
		else:
			prefix = ''
		with open(prefix+'_'+'DTASelect_forwardpeptides','w') as f:
			for peptide in forward_peptide_matches:
				f.write(peptide+'\n')

	print('Finished')

def get_peptides(chunks_with_peptide_matches, peptide_only = False, unique=False):

	''' 
	input: chunks of DTASelect-filter.txt data, one protein ID per chunk
	output: list of peptide sequences (by default, no filtering by unique peptide matches)
	'''

	peptide_matches = []
	for dta_chunk in chunks_with_peptide_matches:
		PSMs = dta_chunk[1:]
		sequence_matches = [PSM.split('\t')[-1] for PSM in PSMs]
		for seq in sequence_matches:
			peptide_matches.append(seq)

	print(len(peptide_matches), 'peptide sequences parsed')

	if peptide_only:
		peptide_matches = [match.split('.')[1] for match in peptide_matches]

	if unique:
		return list(set(peptide_matches))
	else:
		return peptide_matches

def filter_chunks_to_forward_only(chunks_with_peptide_matches, reverse_label='Reverse_'):
	
	'''
	input: protein match chunks from DTASelect-filter.txt
	output: protein match chunks only for 'forward' (true) proteins 
	'''

	forward_chunks = []
	for chunk in chunks_with_peptide_matches:
		if not chunk[0].split('\t')[0].startswith(reverse_label):
			forward_chunks.append(chunk)

	return forward_chunks

def retrieve_loci_noProtDB_old(DTASelect_data):

	'''
	Takes DTASelect data block (for DTASelect-filter.txt file from data run with ProtDB OFF)
	and returns a tuple with length 2

	tuple[0] = list of 'forward' protein loci (ProtDB_IDs)
	tuple[1] = list of 'reverse' protein loci (ProtDB IDs)
	'''
	loci = [line for line in DTASelect_data[2:] if not line.startswith('\t') and not line.startswith('*\t')]

	loci_IDs = [locus_line.split('\t')[0] for locus_line in loci]

	forward_loci = [int(ID) for ID in loci_IDs if not ID.startswith('Reverse_')]
	reverse_loci = [int(ID.replace('Reverse_','')) for ID in loci_IDs if ID.startswith('Reverse_')]

	return (forward_loci, reverse_loci)

def retrieve_loci_noProtDB(DTASelect_data):

	'''
	Takes DTASelect data block (for DTASelect-filter.txt file from data run with ProtDB OFF)
	and returns chunks consisting of an identified locus and peptides from that locus (Forward and Reverse loci)

	One chunk (as an example):
	['61534251\t2\t15\t0.3%\t4485\t511792\t5.8\tU\t8.204104E-6\t0.0069316626\tno description', 
	'\t102214_SC_HEK293_25ug_HCD_FTMS_MS2_12.22000.22000.2\t1.7772602\t0.3232622\t99.2\t1588.8977\t1588.888\t88.0\t1\t10.185568\t11.1\t6\tR.HLSKLFDSLCKLK.F', 
	'\t102214_SC_HEK293_25ug_HCD_FTMS_MS2_12.19162.19162.3\t2.3280215\t0.2957453\t99.8\t1589.9015\t1588.888\t88.0\t1\t8.6662245\t5.6\t9\tR.HLSKLFDSLCKLK.F']

	Generator function
	'''

	loci = []
	current_locus = []
	for line in DTASelect_data[2:]: #skips column names on first 2 lines of data chunk
		if not line.startswith('\t') and not line.startswith('*\t') and current_locus:
			yield current_locus
			current_locus = []
			current_locus.append(line)
		else:
			current_locus.append(line)
	else:
		yield current_locus

def read_DTASelect_filter_file(filepath):

	'''
	Reads in DTASelect-filter.txt file from filepath and returns a 3-item tuple

	tuple[0] = header information (list)
	tuple[1] = file contents list (loci, peptide matches)
	tuple[2] = summary information from footer (dict)
		tuple[2] dict:
		{}
	'''

	with open(filepath) as f:
		DTASelect_file_lines = [line.rstrip('\n') for line in f.readlines()]

	for num, line in enumerate(DTASelect_file_lines):
		if line.startswith('Locus'):
			data_start_num = num
			break

	header = DTASelect_file_lines[:data_start_num]
	header.remove('')
	results = DTASelect_file_lines[data_start_num:-9]
	footer = parse_footer(DTASelect_file_lines[-9:])

	return (header, results, footer)

def parse_footer(footer_lines):

	'''
	Parses list representing this text:
	'
		Proteins	Peptide IDs	Spectra
	Unfiltered	257412	84097	619975
	Filtered	2466	15203	39868
	Forward matches	2423	15117	39745
	Decoy matches	43	86	123
	Forward FDR	1.77	0.57	0.31

	Classification	Nonredundant Proteins	Redundant Proteins
	Unclassified	0	0
	'
	--> ignores last 3 lines for now
	'''

	footer_dict = {}
	table1_keys = footer_lines[0].strip().split('\t')

	for line in footer_lines[1:6]:
		line_split = line.split('\t')
		field = line_split[0]

		line_dict = {}
		for num, column_name in enumerate(table1_keys):
			if field in ('Unfiltered', 'Filtered', 'Forward matches', 'Decoy matches'):
				typecast = int
			else:
				typecast = float
			line_dict[column_name] = typecast(line_split[num+1])

		footer_dict[field] = line_dict

	return footer_dict

if __name__ == '__main__':
	main()
