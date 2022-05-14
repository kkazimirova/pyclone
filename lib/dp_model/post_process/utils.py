import os
from collections import defaultdict

import bz2
import csv

def load_cellular_frequencies_trace(file_name, burnin, thin):
    return _load_trace(file_name, burnin, thin, float)


def load_agregate_cellular_frequencies(analyse_dir, burnin, thin):
    trace = defaultdict(list)

    for dir in os.listdir(analyse_dir):
        trace_file = os.path.join(analyse_dir, dir, 'cellular_frequencies.tsv.bz2')

        fh = bz2.BZ2File(trace_file)
        reader = csv.DictReader(fh, delimiter='\t')

        for i, row in enumerate(reader):
            if i < burnin:
                continue

            if i % thin == 0:
                for mutation in row:
                    trace[mutation].append(float(row[mutation]))  # dict {id_mutation:[values]}

        fh.close()

    return trace



def _load_trace(trace_file, burnin, thin, cast_func):
    trace = defaultdict(list)
    
    fh = bz2.BZ2File(trace_file)
    reader = csv.DictReader(fh, delimiter='\t')


    for i, row in enumerate(reader):
        if i < burnin:
            continue
        
        if i % thin == 0:
            for mutation in row:
                trace[mutation].append(cast_func(row[mutation]))   # dict {id_mutation:[values]}
        
    fh.close()
    
    return trace
