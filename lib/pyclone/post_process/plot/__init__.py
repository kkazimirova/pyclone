import math
import os.path
from collections import OrderedDict

from pyclone.post_process.utils import load_cellular_frequencies_trace, load_cluster_labels_trace

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
    iters = int(math.ceil(mutatations_count / max_size))

    sorter = DensitySorter(trace)
    sorted_clusters = sorter.sort_clusters()
    sorted_clusters = sorted_clusters[0]
    sorted_clusters = list(reversed(sorted_clusters))

    for i in range(iters):
        iter_trace = OrderedDict()
        for j in range(max_size):
            mut = sorted_clusters[(i+1) * max_size - j - 1]
            iter_trace[mut] = trace[mut]

        plot_file = os.path.join(plot_dir, ("plot" + str(i)))
        plotter = CellularFrequencyPlot(iter_trace, cmap=colors)

        print("Plotting " + str(i) + " file of " + str(iters))

        plotter.plot()

        plotter.save(plot_file)


def plot_cellular_frequencies(trace_file, plot_file, burnin, thin):
    trace = load_cellular_frequencies_trace(trace_file, burnin, thin)

    plotter = CellularFrequencyPlot(trace, cmap=colors)

    plotter.plot()

    plotter.save(plot_file)


def my_plot(trace_file, plot_dir, max_size, burnin, thin):
    trace = load_cellular_frequencies_trace(trace_file, burnin, thin)
    mutation_count = len(trace)
    iters = int(math.ceil(mutation_count / max_size))

    sorter = DensitySorter(trace)
    sorted_clusters = sorter.sort_clusters()
    sorted_clusters = sorted_clusters[0]
    sorted_clusters = list(reversed(sorted_clusters))

    for i in range(iters):
        iter_trace = OrderedDict()
        for j in range(max_size):
            mut = sorted_clusters[(i + 1) * max_size - j - 1]
            iter_trace[mut] = trace[mut]

        plot_file = os.path.join(plot_dir, ("plot" + str(i)))
        plotter = MyCellularFrequencyPlot()

        plotter.plot(iter_trace, plot_file)


