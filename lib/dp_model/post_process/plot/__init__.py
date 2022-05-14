import math
import os.path
from collections import OrderedDict

import brewer2mpl

from cellular_frequencies import DensitySorter, MyCellularFrequencyPlot
from densities import PosteriorDensity
from dp_model.post_process.utils import load_cellular_frequencies_trace, load_agregate_cellular_frequencies

bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)

colors = bmap.mpl_colormap


def plot(trace_file, plot_path, max_size, burnin, thin, split=True):
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
            print("Plotting " + str(i + 1) + ". file of " + str(iters) + "...")

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
        print ("Cellular frequencies plotted.")


def agregate_and_plot(analyse_dir, plot_path, max_size, burnin, thin, split=True):
    trace = load_agregate_cellular_frequencies(analyse_dir, burnin, thin)
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

            print ("Plotting " + "plot_" + str(i) + ".pdf...")

            plotter = MyCellularFrequencyPlot()
            plotter.plot(iter_trace, plot_file)

        print ("Plotting done.")

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
        print ("Plotting done.")
