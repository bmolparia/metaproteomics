"""
Merge sqt out files. Combining results
For when, say, searching multiple ptms on the same file
"""

"""
in1 = "/home/gstupp/metaproteomics/test/102214_SC_HEK293_25ug_HCD_FTMS_MS2_10_acetyl.sqt"
in2 = "/home/gstupp/metaproteomics/test/102214_SC_HEK293_25ug_HCD_FTMS_MS2_10_phospho.sqt"
in3 = "/home/gstupp/metaproteomics/test/102214_SC_HEK293_25ug_HCD_FTMS_MS2_10_none.sqt"
"""

#%%
def merge(*args):
    """ recursively merge the sqt files. 
    Merge the first two files in the list, merge that result with the next sqt file, and so on.
    Make sure all inputs are sorted by scan
    """
    #print(len(args))
    if len(args) <= 1:
        raise ValueError
    elif len(args) == 2:
        return merge_two(*args)
    else:
        return merge(merge_two(*args[:2]), *args[2:])

def merge_two(left, right):
    result = []
    while len(left)>0 and len(right)>0:
        new_scan = {'#seq_matching': 0, 'Xcorr': 0, 'deltCN': 0, 'charge': 0, 'high_scan': 0, 'low_scan': 0,
                'lowest_SP': 0, 'matches': [], 'obs_mass': 0, 'process_time': 0, 'server': '', 'total_intensity': 0}
        if left[0]['high_scan'] == right[0]['high_scan']:
            new_scan['matches'] = sorted(left[0]['matches'] + right[0]['matches'], key = lambda x:x['Xcorr'], reverse = True)[:2]
            new_scan['high_scan'] = left[0]['high_scan']
            new_scan['low_scan'] = left[0]['low_scan']
            new_scan['charge'] = left[0]['charge']
            new_scan['matches'][0]['rank_Sp'] = 1
            new_scan['matches'][0]['rank_Xcorr'] = 1
            new_scan['matches'][1]['rank_Sp'] = 2
            new_scan['matches'][1]['rank_Xcorr'] = 2
            new_scan['matches'][0]['deltCN'] = 0
            try:
                new_scan['matches'][1]['deltCN'] = (new_scan['matches'][0]['Xcorr'] - new_scan['matches'][1]['Xcorr']) / new_scan['matches'][0]['Xcorr']
            except ZeroDivisionError:
                new_scan['matches'][1]['deltCN'] = 0
            left.pop(0)
            right.pop(0)
            result.append(new_scan)
        elif left[0]['high_scan'] > right[0]['high_scan']:
            result.append(right.pop(0))
        else:
            result.append(left.pop(0))
    if len(left)>0:
        result.extend(left)
    elif len(right)>0:
        result.extend(right)
    return result
    
def sqt_writer(sqt_chunks, filename):
    spectrum_keys = ['low_scan','high_scan','charge','process_time','server','obs_mass','total_intensity','lowest_SP','#seq_matching']
    match_keys = ['rank_Xcorr','rank_Sp','calc_mass','deltCN','Xcorr','Zscore','matched_ions','expected_ions','seq_matched','valid_status']
    with open(filename, 'w') as f_out:
        for sqt_chunk in sqt_chunks:
            S_line = 'S\t' + '\t'.join([str(sqt_chunk[k]) for k in spectrum_keys]) + '\tU'
            f_out.write(S_line + '\n')
            for match in sqt_chunk['matches']:
                M_line = 'M\t' + '\t'.join([str(match[k]) for k in match_keys])
                f_out.write(M_line + '\n')
                L_lines = '\n'.join(['L\t' + str(L) if not match['reverse'] else 'L\tReverse_' + str(L) for L in match['L']])
                f_out.write(L_lines + '\n')

#%%
import glob
import os
import sys
from file_processing import blazmass_tools
sqt_folder_path = os.getcwd()
sqt_files = glob.glob(os.path.join(sqt_folder_path, "*.sqt"))
# chop off the last "_", merge those
sqt_files = set(['_'.join(x.split('_')[:-1]) for x in sqt_files])
sqt_files.discard('')
for sqt_file in sqt_files:
    sqt_file_paths = glob.glob(sqt_file + "_*.sqt")
    sqt_files = [sorted(list(blazmass_tools.sqt_chunker(sqt_file_path)), key = lambda x:x['low_scan']) for sqt_file_path in sqt_file_paths]
    result = merge(*sqt_files)
    out_file = os.path.join(sqt_folder_path, sqt_file + '.sqt')
    print(out_file)
    sqt_writer(result, out_file)