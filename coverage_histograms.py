import pandas as pd
import numpy as np
import matplotlib as mpl

# to prevent display weirdness when running in Pando:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

HISTOGRAM_BINS = 50
DATA_DIR = './new_models/results/bed_files/'


chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY']
split_index = 11 # ie. chr12
test_chroms = chromosomes[split_index:]

data = pd.read_pickle('./combined_dataset_union_fstitchtfit_with_nascent.pkl')
data = data[(data['chrom'].isin(test_chroms)) & (data['sample'] != 'HCT116')]

positives = data[(data['ovlp_txn'] == 1)]['mean_nr_nascent_reads'].values
plt.clf()
plt.hist(positives[positives <= 0.3], bins=HISTOGRAM_BINS, edgecolor='black', color='green', linewidth=1)
plt.title('Distribution of nascent transcription coverage for peaks \n labeled as "positive". n=%d' % len(positives))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.3))
plt.tight_layout()
plt.savefig("%s/nascent_positives_histogram.png" % DATA_DIR, dpi=300)

negatives = data[(data['ovlp_txn'] == 0)]['mean_nr_nascent_reads'].values
plt.clf()
plt.hist(negatives[negatives <= 0.3], bins=HISTOGRAM_BINS, edgecolor='black', color='green', linewidth=1)
plt.title('Distribution of nascent transcription coverage for peaks \n labeled as "negative". n=%d' % len(negatives))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.3))
plt.tight_layout()
plt.savefig("%s/nascent_negatives_histogram.png" % DATA_DIR, dpi=300)

tp_nascent_coverage = []
fp_nascent_coverage = []
tn_nascent_coverage = []
fn_nascent_coverage = []

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
    tp_nascent_coverage = np.append(tp_nascent_coverage, merged['mean_nr_nascent_reads'].values)


    fp_file = "%s/%s_RNN-hybrid_false_positives.bed" % (DATA_DIR, sample)
    fp_df = pd.read_csv(fp_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, fp_df, on=['chrom', 'start', 'end'])
    fp_nascent_coverage = np.append(fp_nascent_coverage, merged['mean_nr_nascent_reads'].values)


    tn_file = "%s/%s_RNN-hybrid_true_negatives.bed" % (DATA_DIR, sample)
    tn_df = pd.read_csv(tn_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, tn_df, on=['chrom', 'start', 'end'])
    tn_nascent_coverage = np.append(tn_nascent_coverage, merged['mean_nr_nascent_reads'].values)


    fn_file = "%s/%s_RNN-hybrid_false_negatives.bed" % (DATA_DIR, sample)
    fn_df = pd.read_csv(fn_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, fn_df, on=['chrom', 'start', 'end'])
    fn_nascent_coverage = np.append(fn_nascent_coverage, merged['mean_nr_nascent_reads'].values)


tp_nascent_coverage = np.array(tp_nascent_coverage)
fp_nascent_coverage = np.array(fp_nascent_coverage)
tn_nascent_coverage = np.array(tn_nascent_coverage)
fn_nascent_coverage = np.array(fn_nascent_coverage)

plt.clf()
bin_values = plt.hist(tp_nascent_coverage[tp_nascent_coverage < 0.3], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of nascent transcription coverage for all \n true positives. n=%d' % len(tp_nascent_coverage))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.3))
plt.text(0.25, max(bin_values[0]) - max(bin_values[0])/10, "TP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/TP_nascent_histogram.png" % DATA_DIR, dpi=300)

plt.clf()
bin_values = plt.hist(fp_nascent_coverage[fp_nascent_coverage < 0.005], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of nascent transcription coverage for all \n false positives. n=%d' % len(fp_nascent_coverage))
plt.xlabel('Normalized number of reads')
plt.xlim((0,0.005))
plt.ylabel('Count')
plt.text(0.004, max(bin_values[0]) - max(bin_values[0])/10, "FP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/FP_nascent_histogram.png" % DATA_DIR, dpi=300)

plt.clf()
bin_values = plt.hist(tn_nascent_coverage[tn_nascent_coverage < 0.005], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of nascent transcription coverage for all \n true negatives. n=%d' % len(tn_nascent_coverage))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.005))
plt.text(0.004, max(bin_values[0]) - max(bin_values[0])/10, "TN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/TN_nascent_histogram.png" % DATA_DIR, dpi=300)

plt.clf()
bin_values = plt.hist(fn_nascent_coverage[fn_nascent_coverage < 0.3], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of nascent transcription coverage for all \n false negatives. n=%d' % len(fn_nascent_coverage))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.3))
plt.text(0.25, max(bin_values[0]) - max(bin_values[0])/10, "FN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/FN_nascent_histogram.png" % DATA_DIR, dpi=300)



# ATAC-seq

positives = data[(data['ovlp_txn'] == 1)]['mean_nr_reads'].values
plt.clf()
plt.hist(positives[positives < 1.0], bins=HISTOGRAM_BINS, edgecolor='black', color='green', linewidth=1)
plt.title('Distribution of ATAC-seq coverage for peaks \n labeled as "positive". n=%d' % len(positives))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,1.0))
#plt.ylim((0,28000))
plt.tight_layout()
plt.savefig("%s/atac_positives_histogram.png" % DATA_DIR, dpi=300)
del positives

negatives = data[(data['ovlp_txn'] == 0)]['mean_nr_reads'].values
plt.clf()
plt.hist(negatives[negatives < 1.0], bins=HISTOGRAM_BINS, edgecolor='black', color='green', linewidth=1)
plt.title('Distribution of ATAC-seq coverage for peaks \n labeled as "negative". n=%d' % len(negatives))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,1.0))
#plt.ylim((0,28000))
plt.tight_layout()
plt.savefig("%s/atac_negatives_histogram.png" % DATA_DIR, dpi=300)
del negatives

tp_coverage = []
fp_coverage = []
tn_coverage = []
fn_coverage = []

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
    tp_coverage = np.append(tp_coverage, merged['mean_nr_reads'].values)


    fp_file = "%s/%s_RNN-hybrid_false_positives.bed" % (DATA_DIR, sample)
    fp_df = pd.read_csv(fp_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, fp_df, on=['chrom', 'start', 'end'])
    fp_coverage = np.append(fp_coverage, merged['mean_nr_reads'].values)


    tn_file = "%s/%s_RNN-hybrid_true_negatives.bed" % (DATA_DIR, sample)
    tn_df = pd.read_csv(tn_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, tn_df, on=['chrom', 'start', 'end'])
    tn_coverage = np.append(tn_coverage, merged['mean_nr_reads'].values)


    fn_file = "%s/%s_RNN-hybrid_false_negatives.bed" % (DATA_DIR, sample)
    fn_df = pd.read_csv(fn_file, header=None, \
                        sep="\t", na_filter=False, \
                        usecols=[0, 1, 2], \
                        names=['chrom', 'start', 'end'], \
                        dtype={ 'chrom':'str', 'start':'int', \
                                'end  ':'int'})
    merged = pd.merge(data, fn_df, on=['chrom', 'start', 'end'])
    fn_coverage = np.append(fn_coverage, merged['mean_nr_reads'].values)


tp_coverage = np.array(tp_coverage)
fp_coverage = np.array(fp_coverage)
tn_coverage = np.array(tn_coverage)
fn_coverage = np.array(fn_coverage)

plt.clf()
bin_values = plt.hist(tp_coverage[tp_coverage < 0.5], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of ATAC-seq coverage for all \n true positives. n=%d' % len(tp_coverage))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.5))
#plt.ylim((0,6500))
plt.text(0.4, max(bin_values[0]) - max(bin_values[0])/10, "TP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/TP_atac_histogram.png" % DATA_DIR, dpi=300)

plt.clf()
bin_values = plt.hist(fp_coverage[fp_coverage < 0.5], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of ATAC-seq coverage for all \n false positives. n=%d' % len(fp_coverage))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.5))
#plt.ylim((0,6500))
plt.text(0.4, max(bin_values[0]) - max(bin_values[0])/10, "FP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/FP_atac_histogram.png" % DATA_DIR, dpi=300)

plt.clf()
bin_values = plt.hist(tn_coverage[tn_coverage < 0.5], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of ATAC-seq coverage for all \n true negatives. n=%d' % len(tn_coverage))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.5))
#plt.ylim((0,6500))
plt.text(0.4, max(bin_values[0]) - max(bin_values[0])/10, "TN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/TN_atac_histogram.png" % DATA_DIR, dpi=300)

plt.clf()
bin_values = plt.hist(fn_coverage[fn_coverage < 0.5], bins=HISTOGRAM_BINS, edgecolor='black', linewidth=1)
plt.title('Distribution of ATAC-seq coverage for all \n false negatives. n=%d' % len(fn_coverage))
plt.xlabel('Normalized number of reads')
plt.ylabel('Count')
plt.xlim((0,0.5))
#plt.ylim((0,6500))
plt.text(0.4, max(bin_values[0]) - max(bin_values[0])/10, "FN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
plt.tight_layout()
plt.savefig("%s/FN_atac_histogram.png" % DATA_DIR, dpi=300)



