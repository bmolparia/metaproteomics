#%matplotlib inline 

import subprocess
import numpy as np
import glob
from itertools import chain

from bokeh.plotting import show, figure
from bokeh.charts import output_notebook, Bar, output_file


#%% Intensity of PrecursorInt
def log_precursor_int(ms2_file):
    grep = subprocess.Popen(['grep', 'PrecursorInt', ms2_file], stdout=subprocess.PIPE)
    cut = subprocess.Popen('cut -f3'.split(), stdin = grep.stdout, stdout = subprocess.PIPE)
    PrecursorInt = list(map(float,cut.communicate()[0].split()))
    PrecursorInt = [x if x else 1 for x in PrecursorInt]
    LogPrecursorInt = list(map(np.log10, PrecursorInt))
    return LogPrecursorInt

#http://bokeh.pydata.org/en/latest/docs/gallery/histogram.html
def bokeh_precursor_int(LogPrecursorInt):
    p = figure(title = 'Intensity of PrecursorInt')
    output_notebook()
    hist, edges = np.histogram(LogPrecursorInt, density=True, bins=50)
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
         fill_color="#036564", line_color="#033649")
    lpi_median = np.median(LogPrecursorInt)
    p.line([lpi_median, lpi_median], [min(hist),max(hist)], color="#ee3333")
    p.text(text = ["<- median: {0:.2f}".format(lpi_median)], x=lpi_median, y=max(hist)*.9)
    p.xaxis.axis_label = 'log10 intensity'
    p.yaxis.axis_label = 'Count'
    return p

def get_fragment_ions(ms2_file):
    # returns a list of lists
    # each outer list contains list of fragment ion intensities for that scan
    groups = []
    group = []
    numbers = set(list(map(str,range(10))))
    with open(ms2_file) as f:
        for line in f:
            if line[0] in numbers:
                group.append(line.split()[1])
            else:
                if group:
                    groups.append(list(map(np.log10,map(float,group))))
                group = []
        if group:
            groups.append(list(map(np.log10,map(float,group))))
    return groups

#number of fragment ions per scan
def bokeh_num_frag_per_scan(fragment_ions):
    fragment_count = list(map(len,fragment_ions))
    p = figure(title = 'number of fragment ions per scan')
    output_notebook()
    hist, edges = np.histogram(fragment_count, density=True, bins=50)
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
         fill_color="#036564", line_color="#033649")

    fc_median = np.median(fragment_count)
    p.line([fc_median, fc_median], [min(hist),max(hist)], color="#ee3333")
    p.text(text = ["<- median: {0:.2f}".format(fc_median)], x=fc_median, y=max(hist)*.9)
    p.xaxis.axis_label = 'Number of ions'
    p.yaxis.axis_label = 'Count'
    return p

def bokeh_median_frag_int_per_scan(fragment_ions):
    median_int = list(map(np.median,fragment_ions))
    p = figure(title = 'median fragment intensity per scan')
    output_notebook()
    hist, edges = np.histogram(median_int, density=True, bins=50)
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
         fill_color="#036564", line_color="#033649")

    fc_median = np.median(median_int)
    p.line([fc_median, fc_median], [min(hist),max(hist)], color="#ee3333")
    p.text(text = ["<- median: {0:.2f}".format(fc_median)], x=fc_median, y=max(hist)*.9)
    p.xaxis.axis_label = 'Ion Intensity'
    p.yaxis.axis_label = 'Count'
    return p

def bokeh_max_frag_int_per_scan(fragment_ions):
    max_int = list(map(np.max,fragment_ions))
    p = figure(title = 'max fragment intensity per scan')
    output_notebook()
    hist, edges = np.histogram(max_int, density=True, bins=50)
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
         fill_color="#036564", line_color="#033649")

    fc_median = np.median(max_int)
    p.line([fc_median, fc_median], [min(hist),max(hist)], color="#ee3333")
    p.text(text = ["<- median: {0:.2f}".format(fc_median)], x=fc_median, y=max(hist)*.9)
    p.xaxis.axis_label = 'Ion Intensity'
    p.yaxis.axis_label = 'Count'
    return p

def get_lcstep(filename):
    # Parse lc step out of filename
    # failure to parse -> returns -1
    try:
        lcstep = int(filename.split('.')[0].split('_')[-1])
        if lcstep > 100:
            lcstep = int(filename.split('.')[0].split('_')[-2])
        return lcstep
    except ValueError:
        return -1

# Number of scans per ms2 file
def bokeh_num_scans_per_ms2(ms2_list):
    num_scans_per_file = {get_lcstep(ms2_file): int(subprocess.check_output(['grep', '-c', '^S', ms2_file])) for ms2_file in ms2_list}
    x,y = zip(*sorted(num_scans_per_file.items(), key = lambda x:x[0]))
    p = Bar(y, list(map(str,x)), title = 'Number of scans per LC Step', xlabel='LC Step', ylabel='Count')
    return p
#%%
ms2_path = "/mongoc/gstupp/DTASelect/06152015_LC_Optimization/2_fewer_salt_steps_2015_06_16_10_32866"
ms2_list = glob.glob(ms2_path + '/*.ms2')[:3]
LogPrecursorInt = list(chain(*[log_precursor_int(ms2_file) for ms2_file in ms2_list]))
fragment_ions = list(chain(*[get_fragment_ions(ms2_file) for ms2_file in ms2_list]))

show(bokeh_precursor_int(LogPrecursorInt))
show(bokeh_max_frag_int_per_scan(fragment_ions))
show(bokeh_median_frag_int_per_scan(fragment_ions))
show(bokeh_num_frag_per_scan(fragment_ions))
show(bokeh_num_scans_per_ms2(ms2_list))

