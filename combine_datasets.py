#import re
import pickle
import pandas as pd
import numpy as np
import sys


random_labels = False
if len(sys.argv) > 2:
    random_labels = True

DATA_DIR = sys.argv[1]

filenames = [
    '%s/A549-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/GM12878-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/H1-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/HCT116-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/HeLa-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/K562-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/LNCaP-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/MCF7-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
    '%s/THP1-vs-rawnascent_from-fstitchtfit.pk' % DATA_DIR,
]

dframes = dict()
for fname in filenames:
    data = pd.read_pickle(fname)
    atac_signal_cutoff = np.percentile(data.mean_nr_reads.values, 5)
    data = data.drop(data[data['mean_nr_reads'] < atac_signal_cutoff].index)

    dset = fname.split('/')[-1].split('-')[0]
    print(dset)
    # flatten the list of dictionaries into a dataframe and add to our dict
    if random_labels:
        data = data[['chrom', 'start', 'end', 'sequence', 'signal_features', 'mean_nr_reads']]
        labels = [np.random.choice([0,1], p=[0.5, 0.5]) for x in list(range(data.shape[0]))]
        data['bidir_ovlp'] = labels
    data['sample'] = dset
    dframes[dset] = data

data = pd.concat(dframes)

# Try to reduce the data footprint by casting to reasonable datatypes
data['signal_features'] = data['signal_features'].map(lambda x: x.astype('float32'))
data['ovlp_txn'] = data['ovlp_txn'].astype('int')
data['mean_nr_reads'] = data['mean_nr_reads'].astype('float32')
data['gc_ratio'] = data['gc_ratio'].astype('float32')

if random_labels:
    data.to_pickle('/scratch/Users/igtr4848/atac_peak_features/combined_dataset_union_random.pkl')
else:
    data.to_pickle('/scratch/Users/igtr4848/atac_peak_features/combined_dataset_union_fstitchtfit.pkl')

