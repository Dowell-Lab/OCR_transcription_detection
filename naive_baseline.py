import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats
import numpy as np
from sklearn.metrics import confusion_matrix, f1_score, classification_report, roc_curve, auc, accuracy_score
from sklearn.model_selection import train_test_split
import itertools

sns.set(color_codes=True)

DATA_DIR = '.'


data = pd.read_pickle('%s/combined_dataset_union_fstitchtfit_with_nascent.pkl' % DATA_DIR)
data = data[['sample', 'mean_nr_reads', 'mean_nr_nascent_reads', 'ovlp_txn']]
# Change any coverage of zero to a very low number, so that we can still plot the log value
data.loc[(data.mean_nr_reads == 0), 'mean_nr_reads'] = 0.00000001
data.loc[(data.mean_nr_nascent_reads == 0), 'mean_nr_nascent_reads'] = 0.00000001

data['mean_nr_nascent_reads'] = np.log(data.mean_nr_nascent_reads)
data['mean_nr_reads'] = np.log(data.mean_nr_reads)
correlation = np.corrcoef(data['mean_nr_nascent_reads'].values, data['mean_nr_reads'].values)[0,1]
rsq = correlation**2
sns.scatterplot(x="mean_nr_nascent_reads", y="mean_nr_reads", 
                data=data.sort_values('mean_nr_reads'),
                linewidth=0,
                alpha=0.3)
plt.title('Relation between nascent transcription coverage and ATAC-seq coverage\nn=%d  ,  $R^2$=%.6f' % (len(data[(data['mean_nr_reads'] < 10) & (data['mean_nr_nascent_reads'] < 10)]), rsq))
plt.xlabel('log(Normalized mean number of nascent transcription reads)')
plt.ylabel('log(Normalized mean number of ATAC-seq reads)')
plt.xlim((-14, 2))
plt.ylim((-6, 2))
plt.savefig("%s/accessibility-vs-txn_coverage.png" % DATA_DIR, dpi=300)

min_x = min(data.mean_nr_nascent_reads.values)
max_x = max(data.mean_nr_nascent_reads.values)
min_y = min(data.mean_nr_reads.values)
max_y = max(data.mean_nr_reads.values)
palette = itertools.cycle(sns.color_palette())
for cell_type in data['sample'].unique():
    plt.clf()
    ax = sns.scatterplot(x="mean_nr_nascent_reads", y="mean_nr_reads", 
                         data=data[(data['sample'] == cell_type)],
                         linewidth=0,
                         color=next(palette),
                         alpha=0.3)
    ax.set_xlim([min_x, max_x])
    ax.set_ylim([min_y, max_y])
    plt.title('%s cells  ,  n=%d' % (cell_type, len(data[(data['sample'] == cell_type)])), fontsize=20)
    plt.xlabel('log(Mean number of nascent transcription reads)', fontsize=15)
    plt.ylabel('log(Mean number of ATAC-seq reads)', fontsize=15)
    plt.xlim((-14, 2))
    plt.ylim((-6, 2))
    plt.savefig("%s/accessibility-vs-txn_coverage_%s.png" % (DATA_DIR, cell_type), dpi=300)

# Gather the mean number of ATAC-seq reads per OCR where our nascent transcription coverage is
# below the "no transcription" cutoff
negatives = data[(data.ovlp_txn == 0)].mean_nr_reads.values
# Do the same for OCRs where nascent transcription is above the "positive" cutoff
positives = data[(data.ovlp_txn == 1)].mean_nr_reads.values

plt.clf()
sns.distplot(data[(data['mean_nr_reads'] < 0.5)].mean_nr_reads.values, hist=True, kde=True, rug=False, hist_kws={'edgecolor':'black'}, kde_kws={'linewidth': 2, 'alpha':0.5}, bins=100)
plt.title('Distribution of mean number of ATAC-seq reads per OCR')
plt.xlabel('log(ATAC-seq coverage per OCR)')
plt.ylabel('Count')
plt.xlim((-10, 2))
#plt.xlim((-5, 1))
plt.savefig("%s/distribution_atac_coverage.png" % DATA_DIR, dpi=300)

plt.clf()
sns.distplot(data[(data['mean_nr_nascent_reads'] < 0.5)].mean_nr_nascent_reads.values, hist=True, kde=True, rug=False, hist_kws={'edgecolor':'black'}, kde_kws={'linewidth': 2, 'alpha':0.5}, bins=100)
plt.title('Distribution of mean number of nascent transcription reads per OCR')
plt.xlabel('log(Nascent transcription coverage per OCR)')
plt.ylabel('Count')
plt.savefig("%s/distribution_nascent_coverage.png" % DATA_DIR, dpi=300)

plt.clf()
sns.distplot(negatives, hist=True, kde=True, rug=False, hist_kws={'edgecolor':'black'}, kde_kws={'linewidth': 2, 'alpha':0.5}, bins=100)
plt.title('Distribution of mean number of ATAC-seq reads \nfor OCRs without sufficient nascent transcription coverage')
plt.xlabel('log(ATAC-seq coverage per OCR)')
plt.ylabel('Count')
#plt.xlim((-5, 1))
plt.xlim((-10, 2))
plt.savefig("%s/distribution_atac_coverage_negatives.png" % DATA_DIR, dpi=300)


plt.clf()
sns.distplot(positives, hist=True, kde=True, rug=False, hist_kws={'edgecolor':'black'}, kde_kws={'linewidth': 2, 'alpha':0.5}, bins=100)
plt.title('Distribution of mean number of ATAC-seq reads \nfor OCRs with sufficient nascent transcription coverage')
plt.xlabel('log(ATAC-seq coverage per OCR)')
plt.ylabel('Count')
#plt.xlim((-5, 1))
plt.xlim((-10, 2))
plt.savefig("%s/distribution_atac_coverage_positives.png" % DATA_DIR, dpi=300)


X_train, X_test = train_test_split(data, test_size=0.1)
negatives = X_train[(X_train.ovlp_txn == 0)].mean_nr_reads.values
positives = X_train[(X_train.ovlp_txn == 1)].mean_nr_reads.values

print("Training on %d OCRs, testing on %d" % (len(X_train), len(X_test)))
# Use a kernel density estimator with Gaussian kernels to generate a probability distribution
low_kde = stats.gaussian_kde(negatives)
high_kde = stats.gaussian_kde(positives)

y_test = []
y_pred = []
# We can simply calculate the odds ratio for a given OCR, that the mean number of ATAC-seq reads belongs to the "no transcription" vs the "transcription" distributions.
# See:
# https://math.stackexchange.com/questions/825455/probability-that-a-sample-comes-from-one-of-two-distributions
progress = 0
progress_percent = 10
step_size = int(len(X_test)/10)
for ocr in X_test.itertuples():
    progress += 1
    if progress % step_size == 0:
        print("%d%%" % progress_percent)
        progress_percent += 10

    y_test.append(int(ocr.ovlp_txn))

    coverage = ocr.mean_nr_reads
    odds_ratio =  (low_kde.pdf(coverage)/high_kde.pdf(coverage))[0]
    if odds_ratio > 1:
        y_pred.append(0)
    else:
        y_pred.append(1)

np.save('baseline_y_test.npy', y_test)
np.save('baseline_y_pred.npy', y_pred)

print("\n-------------------------------------------------------------------")
print("Baseline Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("-------------------------------------------------------------------")
print("F1-score (macro): %.3f" % f1_score(y_test, y_pred, average='macro'))
print("F1-score (micro): %.3f" % f1_score(y_test, y_pred, average='micro'))
print("F1-score (weighted): %.3f\n" % f1_score(y_test, y_pred, average='weighted'))
print(classification_report(y_test, y_pred))
print("-------------------------------------------------------------------")
print("Accuracy: {}".format(accuracy_score(y_test, y_pred)))

# Plot ROC and report AUC
fpr, tpr, thresholds = roc_curve(y_test, y_pred)
print("AUC: {}".format(auc(fpr, tpr)))

plt.clf()
plt.plot(fpr, tpr);
plt.title('Baseline Reciever Operating Characteristic')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.xlim((0, 1));
plt.ylim((0, 1));
plt.plot(range(2), range(2), '--');
plt.savefig("baseline_ROC.png", dpi=300)

