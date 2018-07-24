"""
Blazmass fails reading a ms2 file with empty scans
Easier (for me) to remove them in python than figure out how blazmass works...

python3 remove_blank_scans path_to_ms2
# warning overwrites input!!!

"""
txt="""S       298     298     1199.01
I       RTime   1.867313
I       BPI     2.366061
I       BPM     1042.795
I       TIC     2.366061
1042.795 2.366061
1043.795 999.3543
10324.23 434.434
S       299     299     1696.97
I       RTime   1.873098
S       300     300     1391.65
I       RTime   1.878788
S       301     301     1231.73
I       RTime   1.884335
I       BPI     5.69242
I       BPM     1760.173
I       TIC     5.69242
1760.173 5.69242
"""
import sys
import string

def get_header(f_handle):
    header = []
    for line in f_handle:
        if line.startswith("H"):
            header.append(line.strip())
        else:
            return '\n'.join(header)
            
def chunker(f_handle):
    lines = []
    for line in f_handle:
        if not line:
            continue
        if line[0] == 'S' and lines:
            yield [line for line in lines if line]
            lines = []
        lines.append(line.strip())
    if lines:
        yield [line for line in lines if line]

scan_is_empty = lambda lines: all(line[0] in string.ascii_uppercase for line in lines)
#%%
with open(sys.argv[1]) as f_handle:
    header = get_header(f_handle)
with open(sys.argv[1]) as f_handle:
    to_write = ['\n'.join(scan) for scan in chunker(f_handle.readlines()) if not scan_is_empty(scan)]


with open(sys.argv[1],'w') as f_handle:
    f_handle.writelines(header)
    f_handle.write('\n')
    f_handle.writelines('\n'.join(to_write))
