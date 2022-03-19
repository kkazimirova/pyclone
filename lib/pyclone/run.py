from __future__ import division

import csv
import math
import os
import random
import shutil

import bz2
import yaml
import pandas

from pyclone.sampler import DirichletProcessSampler, DataPoint
from pyclone.trace import TraceDB
from pyclone.config import load_mutation_from_dict, Mutation, State


def run_dp_model(args):
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


# def cluster_trace(args):
#     from pyclone.post_process.cluster import cluster_pyclone_trace
#
#     pyclone_file = os.path.join(args.trace_dir, 'labels.tsv.bz2')
#
#     print '''Clustering PyClone trace file {in_file} using {method} with a burnin of {burnin} and using every {thin}th sample'''.format(
#         in_file=pyclone_file, method=args.method, burnin=args.burnin, thin=args.thin)
#
#     cluster_pyclone_trace(pyclone_file, args.out_file, args.method, args.burnin, args.thin)



def plot_cellular_frequencies(args):
    import pyclone.post_process.plot as plot

    pyclone_file = os.path.join(args.trace_dir, 'cellular_frequencies.tsv.bz2')


    print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(
        in_file=pyclone_file, burnin=args.burnin, thin=args.thin)

    plot.plot_cellular_frequencies(pyclone_file, args.out_file, args.burnin, args.thin)


def split_and_plot(args):
    import pyclone.post_process.plot as plot

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    pyclone_file = os.path.join(args.trace_dir, 'cellular_frequencies.tsv.bz2')

    # print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(
    #     in_file=pyclone_file, burnin=args.burnin, thin=args.thin)

    plot.split_file_and_plot(pyclone_file, args.out_dir, args.size, args.burnin, args.thin)


# def split_file_and_plot(args):
#     if not os.path.exists(args.out_dir):
#         os.makedirs(args.out_dir)
#
#     size_per_file = 500
#
#     cell_freq_file_path = os.path.join(args.trace_dir, 'cellular_frequencies.tsv.bz2')
#     cell_freq_file = bz2.BZ2File(cell_freq_file_path)
#     df_cell_freq = pandas.read_csv(cell_freq_file, sep='\t')
#
#     columns = len(df_cell_freq.columns)
#
#     iter = int(math.ceil(columns / size_per_file))
#
#     print (columns)
#
#     for i in range(iter):
#         if (i + 1)*size_per_file < columns:
#
#             df_subset = df_cell_freq.iloc[:, i*size_per_file:(i + 1)*size_per_file]
#         else:
#             df_subset = df_cell_freq.iloc[:, i * size_per_file:]
#
#         file_name = 'cell_freq_' + str(i) + '.tsv'
#         path = os.path.join(args.trace_dir, file_name)
#         df_subset.to_csv(path, sep='\t', index=False)
#         print ('created csv')
#         bz2_file = bz2.compress(open(path, 'rb').read())
#         final_file = file_name + '.bz2'
#         final_file_path = os.path.join(args.trace_dir, final_file)
#         bla = open(final_file_path, "wb")
#         bla.write(bz2_file)
#         bla.close()
#         print ('compressed')
#
#         import pyclone.post_process.plot as plot
#         plot_file = os.path.join(args.out_dir, ('plot_' + str(i)))
#         plot.plot_cellular_frequencies(final_file_path, plot_file, args.burnin, args.thin)
#
#     cell_freq_file.close()


# def plot_similarity_matrix(args):
#     import pyclone.post_process.plot as plot
#
#     pyclone_file = os.path.join(args.trace_dir, 'labels.tsv.bz2')
#
#     print '''Plotting similarity matrix from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(
#         in_file=pyclone_file, burnin=args.burnin, thin=args.thin)
#
#     plot.plot_similarity_matrix(pyclone_file, args.out_file, args.burnin, args.thin)


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


def build_random_samples_input_files(args):

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


def random_samples_analyse(args):

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


def random_samples_plot_cf(args):
    import pyclone.post_process.plot as plot

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    for dir in os.listdir(args.in_dir):
        print (dir)

        cellular_file = os.path.join(args.in_dir, dir, 'cellular_frequencies.tsv.bz2')
        out_file = os.path.join(args.out_dir, dir)

        print (cellular_file)
        print (out_file)

        print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(
            in_file=cellular_file, burnin=args.burnin, thin=args.thin)

        plot.plot_cellular_frequencies(cellular_file, out_file, args.burnin, args.thin)

