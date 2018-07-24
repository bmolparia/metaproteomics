#!/usr/bin/env python3

# mass_correct.py
#
# For correcting mass drift in MS2 files
# (if using Garibaldi, need to load Python3 module using 'module load python/3.3.2')
# 
# Sandip Chatterjee
# v1, December 2, 2014
#
# usage: python3 mass_correct.py [-h] [-l label] [-f MS2_input_file.ms2] Xppm
# (corrects all ions in file(s) by X ppm, where X is a float)

import glob
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ppm', help='# of ppm by which to correct masses', type=float)
parser.add_argument('-f', '--ms2file', help='Path to input MS2 file (without agument, will correct all .ms2 files in current directory)')
parser.add_argument('-l', '--label', help='Label for corrected files. For example, File1.ms2 --> File1_corrected.ms2')
parser.add_argument('-p', '--preconly', help='Correct masses for precursor ions only (not fragment ions)', action='store_true')


def main():

	args = parser.parse_args()

	if args.label:
		label = args.label
	else:
		label = '_corrected'

	MS2_files = []
	if args.ms2file:
		MS2_files.append(args.ms2file)
	else:
		print('\nUsing all MS2 files in current directory...')
		MS2_files = glob.glob('*.ms2')
		if not MS2_files:
			print('No MS2 files found. Exiting')
			sys.exit(1)
#	if args.preconly:
#		preconly = True
#	else:
#		preconly = False
	print()
	for MS2_file in MS2_files:
		print('Correcting file {}...'.format(MS2_file))
		correct_MS2_file(MS2_file, args.ppm, args.preconly, label)

	print('Finished')

def get_H_lines(ms2_filename):
	with open(ms2_filename) as f:
		H_lines = [line for line in f if line[0] == 'H']
	return H_lines

def ms2_file_parser(ms2_filename):
    
    ''' returns a list of "chunks" of the MS2 file, with one chunk per scan '''
    ''' (generator function) '''

    current_chunk = []
    with open(ms2_filename) as f:
        for line in f:
            if line[0] != 'H':
                if line[:2] == 'S\t' and current_chunk:
                    yield current_chunk
                    # do something with scan data block (current_chunk) here!

                    current_chunk = []  ## reset current block
                    current_chunk.append(line)
                else:
                    current_chunk.append(line)
        else:
            yield current_chunk

def correct_MS2_file(MS2_file, ppm, preconly=False, label='_corrected'):

	H_lines = get_H_lines(MS2_file)	
	ms2_parser = ms2_file_parser(MS2_file)

	with open(MS2_file.replace('.ms2',label+'.ms2'), 'w') as nf:
		nf.write('H\tMassCorrection\t{0:.2f}ppm\n'.format(ppm))		
		for line in H_lines:
			nf.write(line)
		for scan in ms2_parser:
			corrected_scan = correct_scan(scan, ppm, preconly)
			for line in corrected_scan:
				nf.write(line)

def read_file(file):
    with open(file, 'r') as fin:
        for line in fin:
            yield line.rstrip()


def correct_MS1_file(MS1_file, ppm, label='_corrected'):
    with open(MS1_file.replace('.ms1', label+'.ms1'), 'w') as fout:    
        fout.write('H\tMassCorrection\t{0:.2f}ppm\n'.format(ppm))		        
        for line in read_file(MS1_file):
            if line.startswith(('H', 'S', 'I')):
                fout.write(line+'\n')
            else:
                new_line = line.split(' ')
                new_line[0] = '%1.4f' % get_new_mass(float(new_line[0]), ppm)
                fout.write(' '.join(new_line)+'\n')


def correct_scan(scan, ppm, preconly):

	Z_line = [line for line in scan if line[0] == 'Z'][0]  ## should have one and only one Z line per scan
	S_line = [line for line in scan if line[0] == 'S'][0]  ## should have one and only one S line per scan

	charge, mass = int(Z_line.split('\t')[1]), float(Z_line.split('\t')[2])
	mz = float(S_line.split('\t')[3])

	corrected_mass = get_new_mass(mass, ppm)
	corrected_mz = get_new_mass(mz, ppm)

	corrected_scan = []
	for line in scan:
		if line[0] == 'Z':
			old_Z_line_split = Z_line.split('\t')
			new_Z_line = '\t'.join([old_Z_line_split[0], old_Z_line_split[1], '%.5f' % corrected_mass + '\n'])
			# For example:
			# Z       2       809.59118
			corrected_scan.append(new_Z_line)
		elif line[0] == 'S':
			old_S_line_split = S_line.split('\t')
			new_S_line = '\t'.join([old_S_line_split[0], old_S_line_split[1], old_S_line_split[2], '%.5f' % corrected_mz + '\n'])
			# For example:
			# S       000002  000002  405.29920
			corrected_scan.append(new_S_line)
		elif line[0] == 'I':
			# For example:
			# I       InstrumentType  FTMS
			# I       ActivationType  HCD
			corrected_scan.append(line)
		else:
			# anything else at this point should just be a fragment ion
			if preconly:
				corrected_scan.append(line)
			else:
				# logic to correct fragment
				fragment_mz = float(line.split()[0])
				corrected_fragment_mz = '%.4f' % get_new_mass(fragment_mz, ppm)

				rest_of_fragment_line = ' '.join(line.split()[1:])
				corrected_scan.append(corrected_fragment_mz+' '+rest_of_fragment_line+'\n')

	return corrected_scan

def get_new_mass(obs_mass, ppm):

	'''
	solves for (exact) in the equation:

	ppm = (observed - exact)/exact * 1000000
	'''

	return obs_mass/(1+ppm/1000000)


if __name__ == '__main__':
	main()