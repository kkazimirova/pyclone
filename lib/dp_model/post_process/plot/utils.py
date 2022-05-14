import csv
csv.field_size_limit(10000000)

def load_trace(trace_file, burnin, thin, cast_func):
    trace = {}
    
    reader = csv.DictReader(open(trace_file), delimiter='\t')
    
    for row in reader:
        gene = row['gene']
        
        trace[gene] = [cast_func(x) for x in row['trace'].split(',')][burnin::thin]
    
    return trace