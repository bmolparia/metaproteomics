#!/usr/bin/env python3
# FASTA flattener
# Prints 'flattened' FASTA file to STDOUT
# as in, header, one line seq

import sys

def parse(fasta_file_handle):

    defline, sequence = '', []
    for line in fasta_file_handle:
        if line[0] == '>':
            if defline:
                yield {'defline': defline, 'seq': ''.join(sequence)}
            defline, sequence = line[1:].rstrip(), []
        else:
            sequence.append(line.rstrip())

    if defline:
        yield {'defline': defline, 'seq': ''.join(sequence)}

if __name__ == '__main__':
    for record in parse(sys.stdin):
        print('>' + record['defline'])
        print(record['seq'])
