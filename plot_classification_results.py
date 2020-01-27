import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

my_dpi = 300

data = pd.read_csv("results_using_Adam/results_table.tsv", sep="\t")

sns.set(font_scale=2, style="whitegrid")
g = sns.FacetGrid(data, hue="Cell Type", height=8, hue_kws={"marker": ["o", "p", "s", "P", "d", "^", "X", "*"]} )
sns.violinplot(x="Test OCRs", y="Weighted F1-score", data=data, color="y")
g.map(sns.stripplot, "Test OCRs", "Weighted F1-score", color='black', edgecolor='gray', jitter=0.1, size=15)
plt.savefig("final_results_f1.png", dpi=my_dpi)

plt.clf()
g = sns.FacetGrid(data, hue="Cell Type", height=8, hue_kws={"marker": ["o", "p", "s", "P", "d", "^", "X", "*"]} )
sns.violinplot(x="Test OCRs", y="ROC AUC", data=data, color="c")
g.map(sns.stripplot, "Test OCRs", "ROC AUC", color='black', edgecolor='gray', jitter=0.1, size=15)
plt.savefig("final_results_auc.png", dpi=my_dpi)

plt.clf()
g = sns.FacetGrid(data, hue="Cell Type", height=8, hue_kws={"marker": ["o", "p", "s", "P", "d", "^", "X", "*"]} )
sns.violinplot(x="Test OCRs", y="Training time (minutes)", data=data, color="g")
g.map(sns.stripplot, "Test OCRs", "Training time (minutes)", color='black', edgecolor='gray', jitter=0.1, size=15)
g.add_legend();
plt.savefig("final_results_time.png", dpi=my_dpi)

