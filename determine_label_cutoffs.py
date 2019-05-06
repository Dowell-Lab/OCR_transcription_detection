import pandas as pd
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

cell_types = [
    'gm12878',
    'h1',
    'hct116',
    'k562',
    'cd4pos',
    'jurkat',
    'lncap',
    ]


for cell_type in cell_types:
    data = pd.read_pickle('/scratch/Users/igtr4848/atac_peak_features/%s-vs-rawnascent_from-raw-nascent-reads.pk' % cell_type)
    nascent_cutoff = 0.001
    nascent_cutoff2 = 0.01
    atac_cutoff = np.percentile(data.mean_nr_reads.values, 5)

    plt.clf()
    data2 = data[(data.mean_nr_nascent_reads < 0.02)]
    plt.hist(data2.mean_nr_nascent_reads.values, color='c', bins=50, edgecolor='k')
    plt.axvline(nascent_cutoff, color='k', linestyle='dashed', linewidth=2)
    plt.axvline(nascent_cutoff2, color='k', linestyle='dashed', linewidth=2)
    plt.xlabel('Mean # reads')
    plt.ylabel('Count')
    plt.xlim((0,0.02))
    plt.title('# nascent reads overlaping ATAC-seq peaks for %s cells \n (n = %d)' % (cell_type, len(data)))
    plt.tight_layout()
    plt.savefig("nascent_transcr_coverage_for_%s_peaks.png" % cell_type, dpi=600)

    plt.clf()
    data2 = data[(data.mean_nr_reads < 0.3)]
    plt.hist(data2.mean_nr_reads.values, color='c', bins=50, edgecolor='k')
    plt.xlabel('Mean # reads')
    plt.ylabel('Count')
    plt.xlim((0,0.3))
    plt.title('# ATAC-seq reads overlaping ATAC-seq peaks for %s cells \n (n = %d)' % (cell_type, len(data)))
    plt.tight_layout()
    plt.savefig("atac_transcr_coverage_for_%s_peaks.png" % cell_type, dpi=600)

    del data
