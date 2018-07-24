#! /usr/bin/python

##	fasta_get.py
##	Sandip Chatterjee
##	v3, January 27, 2014
##
##	Usage:
##		for new download:
##		>>python fasta-get.py
##		
##		for resuming directory download (from previous run):
##		>>python fasta-get.py [FILENAME].log
##
##	for obtaining FASTA protein sequence files from the NCBI RefSeq complete database, located at ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/
##	1) downloads FASTA.faa.gz files
##	2) extracts .gz files
##	3) concatenates all files into one large file and removes .gz files (order of files concatenated is not necessarily maintained)

import gzip
from ftplib import FTP
import os
import time
import sys
import shutil

def main():

	startTime = time.time()

	if len(sys.argv) == 1:
		resumeFlag = False
	elif len(sys.argv) == 2:
		resumeFlag = True
		if ".log" not in sys.argv[1]:
			print 'Proper usage: >>fasta_get.py [*.log] ## resume from last downloaded file in this log file'
			sys.exit()
	else:
		print "Wrong number of arguments"
		print 'Proper usage: >>fasta-get.py ## run without resuming'
		print 'Proper usage: >>fasta-get.py *.log ## resume from last downloaded file in this log file'
		sys.exit()

	########################################################################
	####### Change this variable to change RefSeq download directory #######
	########################################################################

	refseq_paths = []
	refseq_paths = [
					'viral',
					'fungi',
					'mitochondrion',
					'plasmid',
					'protozoa',
					'microbial', ##soon to be replaced by the following 2 terms
					# 'archaea',	##available May 2014
					# 'bacteria',	##available May 2014
					]

	refseq_paths = ['complete']	## REMOVE IF NOT DOWNLOADING COMPLETE REFSEQ

	########################################################################
	##### ^ Change this variable to change RefSeq download directory ^ #####
	########################################################################

	refseq_paths_str = '_'.join(refseq_paths)
	refseq_paths = ['/'+x+'/' for x in refseq_paths]

	file_paths = []

	ftp = FTP('ftp.ncbi.nlm.nih.gov')
	ftp.login()

	for refseq_path in refseq_paths:
		
		print 'Accessing NCBI RefSeq'+refseq_path+' at ftp://ftp.ncbi.nlm.nih.gov/refseq/release'+refseq_path
		
		ftp.cwd('/refseq/release'+refseq_path)

		directory_listing = []
		ftp.retrlines("LIST",directory_listing.append)		##	retrieve directory, append to list dirListing

		for directory_line in directory_listing:
			file_paths.append('/refseq/release'+refseq_path+directory_line.split(None, 8)[-1])

	file_paths = [file_path for file_path in file_paths if '.faa.gz' in file_path]

	fasta_master_file = time.strftime("%m%d%y_%H%M%S")+'_refseq_'+refseq_paths_str+'.fasta'

	with open(fasta_master_file, 'a') as extracted_file, open(fasta_master_file+'.log', 'a') as extracted_file_log:
		
		num_files_processed = 0
		
		for count, file_path in enumerate(file_paths):
			
			file_name = file_path.split('/')[-1]	## temporary gzip file

			with open(file_name, 'wb') as temp_gz_file:
				ftp.retrbinary('RETR '+file_path, temp_gz_file.write)

			gzipFile = gzip.open(file_name, 'rb')
			extracted_file.write(gzipFile.read())
			gzipFile.close()

			os.remove(file_name)	##	delete temporary gzip file

			num_files_processed = count+1

			extracted_file_log.write(file_name+" - Item #"+str(num_files_processed)+' of '+str(len(file_paths))+' downloaded, extracted, and appended to master file\n')

		end_status = "Downloaded, extracted, and concatenated "+str(num_files_processed)+" files in "+str((time.time()-startTime)/60)+" minutes\n"
		extracted_file_log.write(end_status)
		extracted_file_log.write("(DATE_CURRENTTIME_refseq_*.fasta)")

	print end_status
	print "Saved compiled data to "+fasta_master_file
	print "(DATE_CURRENTTIME_refseq_*.fasta)"

	ftp.quit()

if __name__ == '__main__':
	main()

		# if resumeFlag:

		# 	##	read in file names from given file log
		# 	logfileNameList = []
		# 	fileLog = open(sys.argv[1], 'rb')
		# 	for line in fileLog:
		# 		logfileNameList.append(line.split()[0])
		# 	logfileNameList = filter(lambda x: '.protein.faa.gz' in x,logfileNameList)	##	only keep lines that contain valid filenames of interest...
		# 	fileLog.close()

		# 	##	concatenate old fasta file and new (empty) fasta file	
		# 	oldFastaFile = open(sys.argv[1].strip('.log'),'rb')
		# 	shutil.copyfileobj(oldFastaFile,extractedFile)
		# 	oldFastaFile.close()

		# 	##	update log file
		# 	fileLog = open(sys.argv[1], 'rb')
		# 	with open(newFileName+'.log', 'wb') as extractedFileLog:
		# 		shutil.copyfileobj(fileLog,extractedFileLog)

		# 	fileLog.close()

		# 	##	remove previously downloaded filenames from new filename (removes order from list due to set operations)
		# 	fileNameList = list(set(fileNameList)-set(logfileNameList))