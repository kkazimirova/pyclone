from __future__ import division

import csv
import os
import random
import shutil
import yaml

from pyclone.sampler import DirichletProcessSampler, DataPoint
from pyclone.trace import TraceDB
from pyclone.config import load_mutation_from_dict, Mutation, State


def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data = load_pyclone_data(args.in_file)

    trace_db = TraceDB(args.out_dir, data.keys())

    try:
        sampler = DirichletProcessSampler(args.tumour_content,
                                          alpha=args.concentration,
                                          alpha_shape=args.concentration_prior_shape,
                                          alpha_rate=args.concentration_prior_rate)
    except:
        trace_db.close()

        shutil.rmtree(args.out_dir)

        raise

    sampler.sample(data.values(), trace_db, num_iters=args.num_iters, seed=args.seed)

    trace_db.close()


def load_pyclone_data(file_name):
    '''
    Load data from PyClone formatted input file.
    '''
    data = {}

    fh = open(file_name)

    config = yaml.load(fh)

    fh.close()

    error_rate = config['error_rate']

    for mutation_dict in config['mutations']:
        mutation = load_mutation_from_dict(mutation_dict)

        data[mutation.id] = DataPoint(mutation.ref_counts,
                                      mutation.var_counts,
                                      mutation.cn_n,
                                      mutation.cn_r,
                                      mutation.cn_v,
                                      mutation.get_mu_n(error_rate),
                                      mutation.get_mu_r(error_rate),
                                      mutation.get_mu_v(error_rate),
                                      mutation.prior_weights)

    return data


def cluster_trace(args):
    from pyclone.post_process.cluster import cluster_pyclone_trace

    pyclone_file = os.path.join(args.trace_dir, 'labels.tsv.bz2')

    print '''Clustering PyClone trace file {in_file} using {method} with a burnin of {burnin} and using every {thin}th sample'''.format(
        in_file=pyclone_file, method=args.method, burnin=args.burnin, thin=args.thin)

    cluster_pyclone_trace(pyclone_file, args.out_file, args.method, args.burnin, args.thin)


def plot_cellular_frequencies(args):
    import pyclone.post_process.plot as plot

    pyclone_file = os.path.join(args.trace_dir, 'cellular_frequencies.tsv.bz2')

    print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(
        in_file=pyclone_file, burnin=args.burnin, thin=args.thin)

    plot.plot_cellular_frequencies(pyclone_file, args.out_file, args.burnin, args.thin)


def plot_similarity_matrix(args):
    import pyclone.post_process.plot as plot

    pyclone_file = os.path.join(args.trace_dir, 'labels.tsv.bz2')

    print '''Plotting similarity matrix from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(
        in_file=pyclone_file, burnin=args.burnin, thin=args.thin)

    plot.plot_similarity_matrix(pyclone_file, args.out_file, args.burnin, args.thin)


def build_input_file(args):
    config = {}

    config['error_rate'] = args.error_rate

    reader = csv.DictReader(open(args.in_file), delimiter='\t')

    config['mutations'] = []

    for row in reader:
        mutation_id = row['mutation_id']

        ref_counts = int(row['ref_counts'])

        var_counts = int(row['var_counts'])

        mutation = Mutation(mutation_id, ref_counts, var_counts)

        cn_n = int(row['cn_n'])

        cn_v = int(row['cn_v'])

        states = _get_states(cn_n, cn_v, args.cn_r, args.g_v)

        for state in states:
            mutation.add_state(state)

        config['mutations'].append(mutation.to_dict())


    fh = open(args.out_file, 'w')

    yaml.dump(config, fh)

    fh.close()


def _get_states(cn_n, cn_v, cn_r_method, g_v_method):
    states = []

    g_v = []

    if g_v_method == 'single':
        g_v.append("A" * (cn_v - 1) + "B")

    elif g_v_method == 'all':
        g_v.append("B" * cn_v)

    elif g_v_method == 'vague':
        for num_var_alleles in range(1, cn_v + 1):
            g_v.append("A" * (cn_v - num_var_alleles) + "B" * num_var_alleles)

    g_n = ["A" * cn_n for _ in g_v]

    if cn_r_method == 'normal':
        g_r = ["A" * cn_n for _ in g_v]

    elif cn_r_method == "variant":
        g_r = ["A" * cn_v for _ in g_v]

    elif cn_r_method == "vague":
        if cn_n == cn_v:
            g_r = ["A" * cn_n for _ in g_v]
        else:
            g_r = ["A" * cn_n for _ in g_v] + ["A" * cn_v for _ in g_v]

            g_n = g_n + g_n

            g_v = g_v + g_v

    prior_weight = [1 for _ in g_v]

    for n, r, v, w in zip(g_n, g_r, g_v, prior_weight):
        states.append(State(n, r, v, w))

    return states


def list_to_csv(l):
    return ",".join([str(x) for x in l])


def klarinka(args):

    out_dir_path = args.out_dir
    if not os.path.exists(out_dir_path):
        os.makedirs(out_dir_path)

    config = {}

    reader = csv.DictReader(open(args.in_file), delimiter='\t')

    config['mutations'] = []

    for row in reader:
        mutation_id = row['mutation_id']

        ref_counts = int(row['ref_counts'])

        var_counts = int(row['var_counts'])

        mutation = Mutation(mutation_id, ref_counts, var_counts)

        cn_n = int(row['cn_n'])

        cn_v = int(row['cn_v'])

        states = _get_states(cn_n, cn_v, args.cn_r, args.g_v)

        for state in states:
            mutation.add_state(state)

        config['mutations'].append(mutation.to_dict())

    for i in range(int(args.count)):
        sample = {}
        sample['error_rate'] = args.error_rate
        sample['mutations'] = random.sample(config['mutations'], int(args.size))

        file_name = '/sample_' + str(i)
        file_path = out_dir_path + file_name

        fh = open(file_path, 'w')

        yaml.dump(sample, fh)

        fh.close()


def klarinka_analyse(args):

    for file_name in os.listdir(args.in_dir):
        print (file_name)

        file_path = args.in_dir + '/' + file_name
        out_dir_path = args.out_dir + '/' + file_name

        print (file_path)
        print (out_dir_path)

        data = load_pyclone_data(file_path)

        trace_db = TraceDB(out_dir_path, data.keys())

        try:
            sampler = DirichletProcessSampler(args.tumour_content,
                                              alpha=args.concentration,
                                              alpha_shape=args.concentration_prior_shape,
                                              alpha_rate=args.concentration_prior_rate)
        except:
            trace_db.close()

            shutil.rmtree(args.out_dir)

            raise

        sampler.sample(data.values(), trace_db, num_iters=args.num_iters, seed=args.seed)

        trace_db.close()
