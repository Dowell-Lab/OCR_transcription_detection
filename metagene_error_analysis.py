import pandas as pd
import numpy as np
from argparse import ArgumentParser
import matplotlib as mpl

# to prevent display weirdness when running in Pando:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

DATA_DIR = './new_models/results/bed_files/'

parser = ArgumentParser(description='This script gets all false positives and false negatives and produces a metagene plot for the ATAC-seq and aggregated nascent transcription signal')
parser.add_argument('-x', '--prefix', dest='prefix', help='Prefix indicating the type of datasets used (histonemarksandtfit, histonemarks, tfitlatest, etc)')
parser.add_argument('-t', '--type', dest='classif_type', help='Classifier type to gather TP and FP from (RNN-hybrid, attr_AdaBoost, sig_SVM, etc)')
args = parser.parse_args()


chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY']
split_index = 11 # ie. chr12
test_chroms = chromosomes[split_index:]

data = pd.read_pickle('./combined_dataset_union_fstitchtfit_with_nascent.pkl')
data = data[(data['chrom'].isin(test_chroms)) & (data['sample'] != 'HCT116')]

mpl.rc('axes', edgecolor='green')
mpl.rc('xtick', color='green')
mpl.rc('ytick', color='green')
positives = data[(data['ovlp_txn'] == 1)]['signal_features'].values
plt.clf()
plt.step(list(range(-500,500)), np.divide(positives.sum(axis=0), len(positives)))
plt.title('ATAC-seq signal meta-peak from all real positives \n n=%d' % len(positives))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total normalized number of reads')
plt.ylim((0,0.22))
plt.tight_layout()
plt.savefig("./ATAC_all_positives_metapeak.png", dpi=300)
del positives

negatives = data[(data['ovlp_txn'] == 0)]['signal_features'].values
plt.clf()
plt.step(list(range(-500,500)), np.divide(negatives.sum(axis=0), len(negatives)))
plt.title('ATAC-seq signal meta-peak from all real negatives \n n=%d' % len(negatives))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total normalized number of reads')
plt.ylim((0,0.22))
plt.tight_layout()
plt.savefig("./ATAC_all_negatives_metapeak.png", dpi=300)
del negatives

tp_atac_signal = []
fp_atac_signal = []
tn_atac_signal = []
fn_atac_signal = []

for sample in data['sample'].unique():
    if sample == 'HCT116':
        # used only for validation
        continue

    tp_file = "%s/%s_RNN-hybrid_true_positives.bed" % (DATA_DIR, sample)
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


    fp_file = "%s/%s_RNN-hybrid_false_positives.bed" % (DATA_DIR, sample)
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


    tn_file = "%s/%s_RNN-hybrid_true_negatives.bed" % (DATA_DIR, sample)
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


    fn_file = "%s/%s_RNN-hybrid_false_negatives.bed" % (DATA_DIR, sample)
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

mpl.rc('axes', edgecolor='black')
mpl.rc('xtick', color='black')
mpl.rc('ytick', color='black')
plt.clf()
plt.step(list(range(-500,500)), np.divide(tp_atac_signal.sum(axis=0), len(tp_atac_signal)))
plt.title('ATAC-seq signal meta-peak from all true positives \n n=%d' % len(tp_atac_signal))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total normalized number of reads')
plt.ylim((0,0.23))
plt.text(400, 0.2, "TP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("./ATAC_TP_metapeak.png", dpi=300)

plt.clf()
plt.step(list(range(-500,500)), np.divide(fp_atac_signal.sum(axis=0), len(fp_atac_signal)))
plt.title('ATAC-seq signal meta-peak from all false positives \n n=%d' % len(fp_atac_signal))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total normalized number of reads')
plt.ylim((0,0.23))
plt.text(400, 0.2, "FP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("./ATAC_FP_metapeak.png", dpi=300)

plt.clf()
plt.step(list(range(-500,500)), np.divide(tn_atac_signal.sum(axis=0), len(tn_atac_signal)))
plt.title('ATAC-seq signal meta-peak from all true negatives \n n=%d' % len(tn_atac_signal))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total normalized number of reads')
plt.ylim((0,0.23))
plt.text(400, 0.2, "TN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("./ATAC_TN_metapeak.png", dpi=300)

plt.clf()
plt.step(list(range(-500,500)), np.divide(fn_atac_signal.sum(axis=0), len(fn_atac_signal)))
plt.title('ATAC-seq signal meta-peak from all false negatives \n n=%d' % len(fn_atac_signal))
plt.xlabel('Nucleotide position from midpoint')
plt.ylabel('Total normalized number of reads')
plt.ylim((0,0.23))
plt.text(400, 0.2, "FN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("./ATAC_FN_metapeak.png", dpi=300)

