import os
## to prevent display weirdness when running in Pando:
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
#plt.ioff()
import numpy as np
import seaborn as sns
import pandas as pd
sns.set(style="whitegrid")

path_to_data = "../rnn_paper_backup/macs2"
path_to_errors = "./new_models/results/bed_files"

samples = [
            'A549',
            'GM12878',
            'H1',
            'HeLa',
            'K562',
            'LNCaP',
            'MCF7',
            'THP1',
          ]

# Proportion of common peaks attributed to TP, FP, etc in entire sample tests

nr_peaks_intersecting_all_samples = float(os.popen("wc -l %s/common_peaks_in_test_set.bed" % path_to_data).read().split()[0])
print(nr_peaks_intersecting_all_samples)

df = pd.DataFrame(columns=['Proportion', 'Cell Type', 'Error Analysis'])
idx = 0

for sample in samples:
    tp_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_true_positives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [tp_overlaps / nr_peaks_intersecting_all_samples, sample, 'True Positives']
    idx += 1

    fp_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_false_positives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [fp_overlaps / nr_peaks_intersecting_all_samples, sample, 'False Positives']
    idx += 1

    tn_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_true_negatives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [tn_overlaps / nr_peaks_intersecting_all_samples, sample, 'True Negatives']
    idx += 1

    fn_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_false_negatives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [fn_overlaps / nr_peaks_intersecting_all_samples, sample, 'False Negatives']
    idx += 1

    print(tp_overlaps, fp_overlaps, tn_overlaps, fn_overlaps)

plt.clf()
my_dpi = 300
fig = plt.figure(figsize=(1300/my_dpi, 800/my_dpi), dpi=my_dpi)
#sns.set(font_scale=2.5)
sns.barplot(x="Proportion", y="Error Analysis", data=df)
plt.title('OCR coordinates overlapped in \n every tested cell type (n=%d)' % nr_peaks_intersecting_all_samples)
plt.tight_layout()
plt.savefig('fig11-common_peak_proportions.png')


# Proportion of TP, FP, etc that are common peaks in entire sample tests

df = pd.DataFrame(columns=['Proportion', 'Cell Type', 'Error Analysis'])
idx = 0

for sample in samples:
    nr_peaks_tp = float(os.popen("wc -l %s/%s_RNN-hybrid_true_positives.bed" % (path_to_errors, sample)).read().split()[0])
    tp_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_true_positives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [tp_overlaps / nr_peaks_tp, sample, 'True Positives']
    idx += 1
    print(sample, nr_peaks_tp, tp_overlaps)

    nr_peaks_fp = float(os.popen("wc -l %s/%s_RNN-hybrid_false_positives.bed" % (path_to_errors, sample)).read().split()[0])
    fp_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_false_positives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [fp_overlaps / nr_peaks_fp, sample, 'False Positives']
    idx += 1
    print(sample, nr_peaks_fp, fp_overlaps)

    nr_peaks_tn = float(os.popen("wc -l %s/%s_RNN-hybrid_true_negatives.bed" % (path_to_errors, sample)).read().split()[0])
    tn_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_true_negatives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [tn_overlaps / nr_peaks_tn, sample, 'True Negatives']
    idx += 1
    print(sample, nr_peaks_tn, tn_overlaps)

    nr_peaks_fn = float(os.popen("wc -l %s/%s_RNN-hybrid_false_negatives.bed" % (path_to_errors, sample)).read().split()[0])
    fn_overlaps = float( os.popen("bedtools intersect -wa -a %s/common_peaks_in_test_set.bed -b %s/%s_RNN-hybrid_false_negatives.bed | wc -l" % (path_to_data, path_to_errors, sample)).read().split()[0])
    df.loc[idx] = [fn_overlaps / nr_peaks_fn, sample, 'False Negatives']
    idx += 1
    print(sample, nr_peaks_fn, fn_overlaps)

plt.clf()
fig = plt.figure(figsize=(3000/my_dpi, 2500/my_dpi), dpi=my_dpi)
#sns.set(font_scale=2.5)
sns.barplot(x="Proportion", y="Error Analysis", data=df)
plt.title('Proportion of peaks in each error analysis category \n attributed to peaks common to all cell types')
plt.tight_layout()
plt.savefig('EA_category_proportion_of_common_peaks.png')

tss = 0
with open("%s/common_peaks_in_test_set.bed" % path_to_data, 'r') as fd:
    for line in fd.readlines():
        chunks = line.split('\t')
        if chunks[-1] != '0\n':
            tss += 1

print("%.2f of OCRs common to all cell types in the test overlap a TSS" % (100*tss/nr_peaks_intersecting_all_samples))

