from pyclone.post_process.utils import load_cellular_frequencies_trace

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == '__main__':
    trace = load_cellular_frequencies_trace('C:\Users\Klara\Documents\PyClone\Pyclone_0.10.0\pyclone\pokusy\ploty\c\cellular_frequencies.tsv.bz2', 0, 1)
    number_of_mutations = len(trace)

    sns.set(style="white")

    plot, axes = plt.subplots(number_of_mutations, 1, figsize=(50, 2 * number_of_mutations))


    plot.tight_layout()


    for i, (key, value) in enumerate(trace.iteritems()):
        data = value
        data = np.array(data)
        subplot = sns.distplot(data, hist=False, kde=True, ax=axes[i], label=key)
        subplot.legend(fontsize=15)
        subplot.set(yticklabels=[])

        l1 = subplot.lines[0]
        x1 = l1.get_xydata()[:, 0]
        y1 = l1.get_xydata()[:, 1]

        subplot.fill_between(x1, y1, color="blue", alpha=0.3)
        # subplot.margins(x=0, y=0.0)

    plot.savefig('C:\Users\Klara\Documents\PyClone\Pyclone_0.10.0\pyclone\pokusy\ploty\lala2.png')
