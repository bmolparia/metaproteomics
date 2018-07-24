#!/usr/bin/env python3

# code borrowed from greg's tax parser

import sys

def main():

    deflines = sys.stdin

    for defline in deflines:
        print(fasta_to_organism_refseq(defline.split('|')[-3]))

def fasta_to_organism_refseq(fasta_defline):
    '''    
    refseq defline are retarded and have no standardized way of noting organism with square brackets
    testers:
    string = "[GSEE] tandem repeats [Invertebrate iridescent virus 30]"
    string = "coat protein [Euphorbia mosaic virus - A [Mexico:Yucatan:2004]]" #pID: 745
    string = "[citrate [pro-3S]-lyase] ligase [Vibrio cholerae]"
    '''
    # If there are 8 pipes (coming from fasta file)
    if fasta_defline.count('|') == 8:
        # Split defline by '|'. Take the 6th
        txt = fasta_defline.split('|')[6].strip()
    elif fasta_defline.count('|') == 4:   
        txt = fasta_defline.split('|')[-1].strip()
    else:
        txt = fasta_defline
    
    try:
        if txt.count('[') != txt.count(']'):
            # logger.warn('Malformed defline: ' + fasta_defline)
             # Just return whatever is between the last '[' and the last ']'   
            txt_flip = txt[::-1]
            organism = txt_flip[txt_flip.find(']')+1:txt_flip.find('[')][::-1]
            return organism
            
    
        brackets = list(parse_brackets(txt))
        # if there are nested brackets
        if max(list(zip(*brackets))[0]) > 0:
            # take the string at level 0
            organism = [x[1] for x in brackets if x[0] == 0]
            return organism[-1]
        else:
            # there are no nested brackets
            # take the last brackets
            return brackets[-1][1]
    except Exception:
        return None

def parse_brackets(string):
    """Generate parenthesized contents in string as pairs (level, contents).
    http://stackoverflow.com/questions/4284991/parsing-nested-parentheses-in-python-grab-content-by-level
    """
    stack = []
    for i, c in enumerate(string):
        if c == '[':
            stack.append(i)
        elif c == ']' and stack:
            start = stack.pop()
            yield (len(stack), string[start + 1: i])

if __name__ == '__main__':
    main()