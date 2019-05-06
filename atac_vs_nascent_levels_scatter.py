import pandas as pd
import numpy as np
import matplotlib as mpl

# to prevent display weirdness when running in Pando:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import seaborn as sns

DATA_DIR = '/scratch/Users/igtr4848/atac_peak_features'

data = pd.read_pickle('%s/combined_dataset_union_rawnascent.pkl' % DATA_DIR)

no_outliars = data[(data['mean_nr_nascent_reads'] < 10) & (data['mean_nr_reads'] < 10)]
del data

correlation = np.corrcoef(no_outliars['mean_nr_nascent_reads'].values, no_outliars['mean_nr_reads'].values)[0,1]
rsq = correlation**2

plt.clf()
sns.scatterplot(x=np.log(no_outliars['mean_nr_nascent_reads'].values), y=np.log(no_outliars['mean_nr_reads'].values))
plt.title('Relation between nascent transcription coverage and ATAC-seq coverage \n n=%d , $R^2$=%.4f' % (len(no_outliars), rsq))
plt.xlabel('log(Normalized mean # nascent reads)')
plt.ylabel('log(Normalized mean # ATAC-seq reads)')
plt.tight_layout()
plt.savefig("%s/accessibility-vs-txn_coverage.png" % DATA_DIR, dpi=300)

