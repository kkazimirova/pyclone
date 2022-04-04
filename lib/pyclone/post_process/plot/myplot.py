import math

from pyclone.post_process.utils import load_cellular_frequencies_trace

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == '__main__':
    sns.set(style="white")

    trace = load_cellular_frequencies_trace('C:\Users\Klara\Documents\PyClone\Pyclone_0.10.0\pyclone\pokusy\ploty\\2.txt', 0, 1)

    number_of_mutations = len(trace)
    plot, axes = plt.subplots(number_of_mutations, 1, figsize=(50, 2 * number_of_mutations))
    plot.tight_layout()

    for i, (key, value) in enumerate(trace.iteritems()):
        data = value
        data = np.array(data)
        subplot = sns.distplot(data, hist=False, kde=True, kde_kws={"bw":0.007}, ax=axes[i], label=key)
        subplot.legend(fontsize=15)
        subplot.set(yticklabels=[])
        axes[i].set_xlim(0, 1)
        axes[i].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

        l1 = subplot.lines[0]
        x1 = l1.get_xydata()[:, 0]
        y1 = l1.get_xydata()[:, 1]

        subplot.fill_between(x1, y1, color="blue", alpha=0.3)


    plot.savefig('C:\Users\Klara\Documents\PyClone\Pyclone_0.10.0\pyclone\pokusy\ploty\porazima2.png')
    print ("la")
    plt.close(plot)






