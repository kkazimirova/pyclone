import math
import os.path
from collections import OrderedDict

from pyclone.post_process.utils import load_cellular_frequencies_trace, load_agregate_cellular_frequencies

from cellular_frequencies import CellularFrequencyPlot, DensitySorter, MyCellularFrequencyPlot
from densities import PosteriorDensity
from similarity_matrix import SimilarityMatrixPlot

import brewer2mpl

bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)

colors = bmap.mpl_colormap


# def plot_similarity_matrix(trace_file, plot_file, burnin, thin):
#     trace = load_cluster_labels_trace(trace_file, burnin, thin)
#
#     if len(trace) < 2:
#         return
#
#     plotter = SimilarityMatrixPlot(trace)
#
#     plotter.plot()
#
#     plotter.save(plot_file)


def split_file_and_plot(trace_file, plot_dir, max_size, burnin, thin):
    trace = load_cellular_frequencies_trace(trace_file, burnin, thin)
    mutatations_count = len(trace)
    iters = int(math.ceil(mutatations_count / float(max_size)))

    sorter = DensitySorter(trace)
    sorted_clusters = sorter.sort_clusters()
    sorted_clusters = sorted_clusters[0]
    sorted_clusters = list(reversed(sorted_clusters))

    for i in range(iters):
        iter_trace = OrderedDict()

        if i == iters - 1:
            for j in range(mutatations_count % max_size):
                mut = sorted_clusters[mutatations_count - j - 1]
                iter_trace[mut] = trace[mut]
        else:
            for j in range(max_size):
                mut = sorted_clusters[(i+1) * max_size - j - 1]
                iter_trace[mut] = trace[mut]

        plot_file = os.path.join(plot_dir, ("plot" + str(i)))
        plotter = CellularFrequencyPlot(iter_trace, cmap=colors)

        print("Plotting " + str(i) + " file of " + str(iters - 1))

        plotter.plot()
        plotter.save(plot_file)


def plot_cellular_frequencies(trace_file, plot_file, burnin, thin):
    trace = load_cellular_frequencies_trace(trace_file, burnin, thin)

    plotter = CellularFrequencyPlot(trace, cmap=colors)

    plotter.plot()
    plotter.save(plot_file)


def my_plot(trace_file, plot_path, max_size, burnin, thin, split=True):
    trace = load_cellular_frequencies_trace(trace_file, burnin, thin)
    mutation_count = len(trace)

    if split:
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)

        iters = int(math.ceil(mutation_count / float(max_size)))

        sorter = DensitySorter(trace)
        sorted_clusters = sorter.sort_clusters()
        sorted_clusters = sorted_clusters[0]
        sorted_clusters = list(reversed(sorted_clusters))

        for i in range(iters):
            iter_trace = OrderedDict()

            if i == iters - 1:
                for j in range(mutation_count % max_size):
                    mut = sorted_clusters[(i * max_size) + j]
                    iter_trace[mut] = trace[mut]
            else:
                for j in range(max_size):
                    mut = sorted_clusters[(i * max_size) + j]
                    iter_trace[mut] = trace[mut]

            plot_file = os.path.join(plot_path, ("plot_" + str(i) + ".pdf"))
            print(plot_file)

            plotter = MyCellularFrequencyPlot()
            plotter.plot(iter_trace, plot_file)

    else:
        sorter = DensitySorter(trace)
        sorted_clusters = sorter.sort_clusters()
        sorted_clusters = sorted_clusters[0]
        sorted_clusters = list(reversed(sorted_clusters))

        iter_trace = OrderedDict()
        for i in range(mutation_count):
            mut = sorted_clusters[i]
            iter_trace[mut] = trace[mut]

        plotter = MyCellularFrequencyPlot()
        plotter.plot(iter_trace, (plot_path + ".pdf"))


def agregate_and_plot(analyse_dir, plot_path, max_size, burnin, thin, split=True):
    print("analyse dir  " + analyse_dir)
    trace = load_agregate_cellular_frequencies(analyse_dir, burnin, thin)
    mutation_count = len(trace)

    for key in trace:
        print (key, len(trace[key]))

    if split:
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)

        iters = int(math.ceil(mutation_count / float(max_size)))

        sorter = DensitySorter(trace)
        sorted_clusters = sorter.sort_clusters()
        sorted_clusters = sorted_clusters[0]
        sorted_clusters = list(reversed(sorted_clusters))

        for i in range(iters):
            iter_trace = OrderedDict()

            if i == iters - 1:
                for j in range(mutation_count % max_size):
                    mut = sorted_clusters[(i * max_size) + j]
                    iter_trace[mut] = trace[mut]
            else:
                for j in range(max_size):
                    mut = sorted_clusters[(i * max_size) + j]
                    iter_trace[mut] = trace[mut]

            plot_file = os.path.join(plot_path, ("plot_" + str(i) + ".pdf"))

            plotter = MyCellularFrequencyPlot()
            plotter.plot(iter_trace, plot_file)

    else:
        sorter = DensitySorter(trace)
        sorted_clusters = sorter.sort_clusters()
        sorted_clusters = sorted_clusters[0]
        sorted_clusters = list(reversed(sorted_clusters))

        iter_trace = OrderedDict()
        for i in range(mutation_count):
            mut = sorted_clusters[i]
            iter_trace[mut] = trace[mut]

        plotter = MyCellularFrequencyPlot()
        plotter.plot(iter_trace, (plot_path + ".pdf"))


