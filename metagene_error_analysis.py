import pandas as pd
import numpy as np
from argparse import ArgumentParser
import matplotlib as mpl

# to prevent display weirdness when running in Pando:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()




DATA_DIR = '/scratch/Users/igtr4848/atac_peak_features'

parser = ArgumentParser(description='This script gets all false positives and false negatives and produces a metagene plot for the ATAC-seq and aggregated nascent transcription signal')
parser.add_argument('-x', '--prefix', dest='prefix', help='Prefix indicating the type of datasets used (histonemarksandtfit, histonemarks, tfitlatest, etc)')
parser.add_argument('-t', '--type', dest='classif_type', help='Classifier type to gather TP and FP from (RNN-hybrid, attr_AdaBoost, sig_SVM, etc)')
args = parser.parse_args()

samples = [
            'gm12878',
            'h1',
            'hct116',
            'k562',
            'cd4pos',
            'jurkat',
            'lncap',
          ]

data = pd.read_pickle('%s/combined_dataset_union_%s.pkl' % (DATA_DIR, args.prefix))

positives = data[(data['mean_nr_nascent_reads'] > 0.01)]['signal_features'].values
plt.clf()
plt.step(list(range(-500,500)), positives.sum(axis=0))
plt.title('ATAC-seq signal meta-peak from all real positives \n n=%d' % len(positives))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total # normalized reads')
plt.tight_layout()
plt.savefig("%s/positives_%s_%s_metapeak.png" % (DATA_DIR, args.prefix, args.classif_type), dpi=300)
del positives

negatives = data[(data['mean_nr_nascent_reads'] < 0.001)]['signal_features'].values
plt.clf()
plt.step(list(range(-500,500)), negatives.sum(axis=0))
plt.title('ATAC-seq signal meta-peak from all real negatives \n n=%d' % len(negatives))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total # normalized reads')
plt.tight_layout()
plt.savefig("%s/negatives_%s_%s_metapeak.png" % (DATA_DIR, args.prefix, args.classif_type), dpi=300)
del negatives

tp_atac_signal = []
fp_atac_signal = []
tn_atac_signal = []
fn_atac_signal = []

for sample in samples:
    tp_file = "%s/%s-vs-rawnascent_RNN-hybrid_true_positives.bed" % (DATA_DIR, sample)
    tp_df = pd.read_csv(tp_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, tp_df, on=['chrom', 'start', 'end'])
    if len(tp_atac_signal) > 0:
        tp_atac_signal = np.append(tp_atac_signal, merged['signal_features'].values)
    else:
        tp_atac_signal = merged['signal_features'].values


    fp_file = "%s/%s-vs-rawnascent_RNN-hybrid_false_positives.bed" % (DATA_DIR, sample)
    fp_df = pd.read_csv(fp_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, fp_df, on=['chrom', 'start', 'end'])
    if len(fp_atac_signal) > 0:
        fp_atac_signal = np.append(fp_atac_signal, merged['signal_features'].values)
    else:
        fp_atac_signal = merged['signal_features'].values


    tn_file = "%s/%s-vs-rawnascent_RNN-hybrid_true_negatives.bed" % (DATA_DIR, sample)
    tn_df = pd.read_csv(tn_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, tn_df, on=['chrom', 'start', 'end'])
    if len(tn_atac_signal) > 0:
        tn_atac_signal = np.append(tn_atac_signal, merged['signal_features'].values)
    else:
        tn_atac_signal = merged['signal_features'].values


    fn_file = "%s/%s-vs-rawnascent_RNN-hybrid_false_negatives.bed" % (DATA_DIR, sample)
    fn_df = pd.read_csv(fn_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, fn_df, on=['chrom', 'start', 'end'])
    if len(fn_atac_signal) > 0:
        fn_atac_signal = np.append(fn_atac_signal, merged['signal_features'].values)
    else:
        fn_atac_signal = merged['signal_features'].values


tp_atac_signal = np.array(tp_atac_signal)
fp_atac_signal = np.array(fp_atac_signal)
tn_atac_signal = np.array(tn_atac_signal)
fn_atac_signal = np.array(fn_atac_signal)

plt.clf()
plt.step(list(range(-500,500)), tp_atac_signal.sum(axis=0))
plt.title('ATAC-seq signal meta-peak from all true positives \n n=%d (%s / %s)' % (len(tp_atac_signal), args.prefix, args.classif_type))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total # normalized reads')
plt.tight_layout()
plt.savefig("%s/TP_%s_%s_metapeak.png" % (DATA_DIR, args.prefix, args.classif_type), dpi=300)

plt.clf()
plt.step(list(range(-500,500)), fp_atac_signal.sum(axis=0))
plt.title('ATAC-seq signal meta-peak from all false positives \n n=%d (%s / %s)' % (len(fp_atac_signal), args.prefix, args.classif_type))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total # normalized reads')
plt.tight_layout()
plt.savefig("%s/FP_%s_%s_metapeak.png" % (DATA_DIR, args.prefix, args.classif_type), dpi=300)

plt.clf()
plt.step(list(range(-500,500)), tn_atac_signal.sum(axis=0))
plt.title('ATAC-seq signal meta-peak from all true negatives \n n=%d (%s / %s)' % (len(tn_atac_signal), args.prefix, args.classif_type))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total # normalized reads')
plt.tight_layout()
plt.savefig("%s/TN_%s_%s_metapeak.png" % (DATA_DIR, args.prefix, args.classif_type), dpi=300)

plt.clf()
plt.step(list(range(-500,500)), fn_atac_signal.sum(axis=0))
plt.title('ATAC-seq signal meta-peak from all false negatives \n n=%d (%s / %s)' % (len(fn_atac_signal), args.prefix, args.classif_type))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total # normalized reads')
plt.tight_layout()
plt.savefig("%s/FN_%s_%s_metapeak.png" % (DATA_DIR, args.prefix, args.classif_type), dpi=300)

