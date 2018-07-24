#! /bin/bash
#"""
#Truncate fasta defline to 30 characters
#"""
awk '{if ($1 !~ /^\w/) {print substr($1, 1, 29)} else {print $1}}'

