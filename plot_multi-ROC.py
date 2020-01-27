import matplotlib as mpl

## to prevent display weirdness when running in Pando:
#mpl.use('Agg')
import matplotlib.pyplot as plt
#plt.ioff()
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import auc, f1_score, precision_recall_curve, roc_curve, roc_auc_score
import itertools

sns.set(color_codes=True)

cell_types = ['A549', 'GM12878', 'H1', 'HeLa', 'K562', 'LNCaP', 'MCF7', 'THP1']

RESULTS_DIR = "./results_using_Adam"

my_dpi = 300

for cell in cell_types:
    plt.clf()
    fig = plt.figure(figsize=(1300/my_dpi, 1200/my_dpi), dpi=my_dpi)
    print('plotting %s TSS' % cell)
    ax = sns.lineplot(x=[0,1], y=[0,1], dashes=[(2, 2), (2, 2)])
    ax.lines[0].set_linestyle("--")
    y_test = np.load('%s/%s_y_test.npy' % (RESULTS_DIR, cell))
    y_prob = np.load('%s/%s_y_prob.npy' % (RESULTS_DIR, cell))
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    roc_auc = roc_auc_score(y_test, y_prob)
    sns.lineplot(x=fpr, y=tpr, label='All')
    tss_y_test = np.load('%s/tss_%s_y_test.npy' % (RESULTS_DIR, cell))
    tss_y_prob = np.load('%s/tss_%s_y_prob.npy' % (RESULTS_DIR, cell))
    tss_fpr, tss_tpr, _ = roc_curve(tss_y_test, tss_y_prob)
    sns.lineplot(x=tss_fpr, y=tss_tpr, label='TSS')
    notss_y_test = np.load('%s/notss_%s_y_test.npy' % (RESULTS_DIR, cell))
    notss_y_prob = np.load('%s/notss_%s_y_prob.npy' % (RESULTS_DIR, cell))
    notss_fpr, notss_tpr, _ = roc_curve(notss_y_test, notss_y_prob)
    sns.lineplot(x=notss_fpr, y=notss_tpr, label='Non-TSS')
    plt.title("Receiver Operating Characteristic \n(%s, AUC=%.2f)" % (cell, roc_auc))
    plt.margins(x=0)
    plt.margins(y=0)
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.savefig('tss_notss_roc_%s.png' % cell)

    plt.clf()
    fig = plt.figure(figsize=(1300/my_dpi, 1200/my_dpi), dpi=my_dpi)
    precision, recall, _ = precision_recall_curve(y_test, y_prob)
    sns.lineplot(x=recall, y=precision, label='All')
    tss_precision, tss_recall, _ = precision_recall_curve(tss_y_test, tss_y_prob)
    sns.lineplot(x=tss_recall, y=tss_precision, label='TSS')
    notss_precision, notss_recall, _ = precision_recall_curve(notss_y_test, notss_y_prob)
    sns.lineplot(x=notss_recall, y=notss_precision, label='Non-TSS')
    y_pred = np.load('%s/%s_y_pred.npy' % (RESULTS_DIR, cell))
    f1 = f1_score(y_test, y_pred, average='weighted')
    plt.title("Precision/Recall \n(%s, F1=%.2f)" % (cell, f1))
    plt.margins(x=0)
    plt.margins(y=0)
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.savefig('tss_notss_pr_%s.png' % cell)
    
plt.clf()
fig = plt.figure(figsize=(1300/my_dpi, 1300/my_dpi), dpi=my_dpi)
palette = itertools.cycle(sns.color_palette())
ax = sns.lineplot(x=[0,1], y=[0,1], dashes=[(2, 2), (2, 2)])
ax.lines[0].set_linestyle("--")
for cell in cell_types:
    print("Processing curve for %s" % cell)
    if cell == 'HeLa':
        # preseve the same colors as the ATAC-vs-nascent coverage plots
        next(palette)
    y_test = np.load('%s/%s_y_test.npy' % (RESULTS_DIR, cell))
    y_prob = np.load('%s/%s_y_prob.npy' % (RESULTS_DIR, cell))
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    roc_auc = roc_auc_score(y_test, y_prob)
    sns.lineplot(x=fpr, y=tpr, label='%s' % cell, color=next(palette))
plt.margins(x=0)
plt.margins(y=0)
plt.xlim((0, 1))
plt.ylim((0, 1))
plt.savefig('fig4-roc_curves_samples.png', dpi=300)

