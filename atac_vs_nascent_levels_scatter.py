import pandas as pd
import numpy as np
import matplotlib as mpl
from scipy.stats import gaussian_kde

# to prevent display weirdness when running in Pando:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import seaborn as sns


DATA_DIR = '.'

data = pd.read_pickle('%s/combined_dataset_union_fstitchtfit_with_nascent.pkl' % DATA_DIR)
data = data[['sample', 'mean_nr_reads', 'mean_nr_nascent_reads', 'ovlp_txn']]

# Change any coverage of zero to a very low number, so that we can still plot the log value
data.loc[(data.mean_nr_reads == 0), 'mean_nr_reads'] = 0.00000001
data.loc[(data.mean_nr_nascent_reads == 0), 'mean_nr_nascent_reads'] = 0.00000001
data['mean_nr_nascent_reads'] = np.log(data.mean_nr_nascent_reads)
data['mean_nr_reads'] = np.log(data.mean_nr_reads)

nascent_reads = data['mean_nr_nascent_reads'].values
atac_reads = data['mean_nr_reads'].values

correlation = np.corrcoef(nascent_reads, atac_reads)[0,1]
rsq = correlation**2

xy = np.vstack([nascent_reads, atac_reads])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
np_nascent_coverage, np_atac_coverage, z = nascent_reads[idx], atac_reads[idx], z[idx]

plt.clf()
fig, ax = plt.subplots()
ax.scatter(x=np_nascent_coverage, y=np_atac_coverage, c=z)
plt.title('Relation between nascent transcription coverage and ATAC-seq coverage \n n=%d , $R^2$=%.6f' % (len(nascent_reads), rsq))
plt.xlabel('log(Normalized mean number of nascent transcription reads)')
plt.ylabel('log(Normalized mean number of ATAC-seq reads)')
plt.xlim((-14.0, 2.0));
plt.ylim((-6.0, 2.0));
plt.tight_layout()
plt.savefig("%s/accessibility-vs-txn_coverage.png" % DATA_DIR, dpi=300)

