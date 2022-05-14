from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


class DensitySorter(object):
    def __init__(self, frequencies_trace):
        self.trace = frequencies_trace
        self._init_clusters()

    def _init_clusters(self):
        self.gene_color = {}

        self.clusters = [self.trace.keys(), ]  # [[id_mutacii]]

    def sort_clusters(self):
        sort_objective = np.mean

        cluster_values = [[] for _ in self.clusters]

        for i, cluster in enumerate(self.clusters):
            for gene in cluster:
                cluster_values[i].extend(self.trace[gene])

        sort_values = []

        for value in cluster_values:
            sort_values.append(sort_objective(value))

        new_clusters = []

        for _, cluster in sorted(zip(sort_values, self.clusters)):
            new_clusters.append(cluster)

        for i, cluster in enumerate(new_clusters):
            sort_values = [sort_objective(self.trace[gene]) for gene in cluster]

            new_clusters[i] = []

            for _, gene in sorted(zip(sort_values, cluster)):
                new_clusters[i].append(gene)

        return new_clusters


class MyCellularFrequencyPlot(object):
    def __init__(self, frequencies_trace=None, clusters=None):
        pass

    def plot(self, trace, out_file, sort_clusters=True, sort_genes=True):
        sns.set(style="white")

        number_of_mutations = len(trace)
        plot, axes = plt.subplots(number_of_mutations, 1, figsize=(30, 2 * number_of_mutations))
        plot.tight_layout()

        for i, (key, value) in enumerate(trace.iteritems()):
            data = value
            data = np.array(data)
            subplot = sns.distplot(data, hist=False, kde=True, kde_kws={"bw": 0.01}, ax=axes[i], label=key)
            subplot.legend(fontsize=15)
            subplot.set(yticklabels=[])
            axes[i].set_xlim(0, 1)
            axes[i].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

            l1 = subplot.lines[0]
            x1 = l1.get_xydata()[:, 0]
            y1 = l1.get_xydata()[:, 1]

            subplot.fill_between(x1, y1, color="blue", alpha=0.3)
            # subplot.margins(x=0, y=0.0)

        plot.savefig(out_file)
        plt.close(plot)
