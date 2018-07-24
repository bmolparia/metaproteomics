#!/usr/bin/env python3
'''
reads FASTA file on STDIN, outputs proteins to STDOUT that don't begin with 'Reverse'
made for the ProtDB FASTA defline format:

# >209||gi|430025802|ref|YP_007195275.1| VP2 [Cebus albifrons polyomavirus 1]
# >82818599||Reverse_gi|446730281|ref|YP_007173651.1| phage protein-Restriction nuclease superfamily [Lactobacillus phage phiAQ113]

# Usage: cat renumbered_071414_indexDB_reversed.fasta | python ~/metaproteomics/remove_reverse_fasta_records.py > renumbered_071414_indexDB_forward.fasta
or 
# cat ~/renumbered_071414_indexDB_reversed.fasta | parallel -j+8 --block 500K --pipe python ~/metaproteomics/remove_reverse_fasta_records.py > renumbered_071414_indexDB_forward.fasta

Greg Stupp

'''
import sys
from Bio import SeqIO
parser = SeqIO.parse(sys.stdin, 'fasta')
for p in parser:
    if not p.description.split('||')[1].startswith('Reverse_'):
        print(p.format('fasta').strip())