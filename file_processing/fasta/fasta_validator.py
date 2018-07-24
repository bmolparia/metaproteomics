#!/usr/bin/env python3
# protein FASTA validator
# Make sure seq only contains valid AAs

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

AA = set("ACDEFGHIKLMNPQRSTVWY")
if __name__ == '__main__':
    for record in parse(sys.stdin):
        if len(set(record['seq']) - AA) > 0:
            continue
        print('>' + record['defline'])
        print(record['seq'])
