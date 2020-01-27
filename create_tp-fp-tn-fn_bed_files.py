import pandas as pd
import numpy as np

DATA_DIR = "./new_models/results/bed_files"
data = pd.read_pickle('./combined_dataset_union_fstitchtfit.pkl')

chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY']
split_index = 11 # ie. chr12
test_chroms = chromosomes[split_index:]

for test_index in range(8):
    train_samples = [
        'A549',
        'GM12878',
        'H1', 
        'HeLa',
        'K562',
        'LNCaP',
        'MCF7',
        'THP1']
    test_sample = train_samples[test_index]

    false_positives = []
    false_negatives = []
    true_positives = []
    true_negatives = []
    y_pred = np.load("./results_using_Adam/%s_y_pred.npy" % test_sample)
    y_test = np.load("./results_using_Adam/%s_y_test.npy" % test_sample)

    i = 0
    for ocr in data[(data['sample'] == test_sample) & (data['chrom'].isin(test_chroms))].itertuples():
        if ocr.ovlp_txn == 0 and y_pred[i] == 1:
            false_positives.append([ocr.chrom, ocr.start, ocr.end])
        elif ocr.ovlp_txn == 1 and y_pred[i] == 0:
            false_negatives.append([ocr.chrom, ocr.start, ocr.end])
        elif ocr.ovlp_txn == 0 and y_pred[i] == 0:
            true_negatives.append([ocr.chrom, ocr.start, ocr.end])
        else:
            true_positives.append([ocr.chrom, ocr.start, ocr.end])
        i += 1

    with open("%s/%s_RNN-hybrid_false_positives.bed" % (DATA_DIR, test_sample), 'w') as fp_file:
        for peak in false_positives:
            fp_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))
    with open("%s/%s_RNN-hybrid_false_negatives.bed" % (DATA_DIR, test_sample), 'w') as fn_file:
        for peak in false_negatives:
            fn_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))
    with open("%s/%s_RNN-hybrid_true_positives.bed" % (DATA_DIR, test_sample), 'w') as tp_file:
        for peak in true_positives:
            tp_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))
    with open("%s/%s_RNN-hybrid_true_negatives.bed" % (DATA_DIR, test_sample), 'w') as tn_file:
        for peak in true_negatives:
            tn_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))

