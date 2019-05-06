import sys
import pandas as pd
import numpy as np
from keras.preprocessing.text import Tokenizer
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import ReduceLROnPlateau, EarlyStopping
from sklearn.metrics import confusion_matrix, f1_score, classification_report, roc_curve, auc, accuracy_score
from sklearn.pipeline import FeatureUnion, Pipeline
from sklearn.model_selection import train_test_split
from FeatureExtraction import ItemSelector, Reshape, Seq2Ind
from rnn_classifier import RNNHybridClassifier

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

NR_GPUS = 1
LR = float(sys.argv[1])
HIDDEN_SIZE = int(sys.argv[2])
TEST_CHROM = sys.argv[3]   # e.g. chr12, chrX, etc
TEST_SAMPLE = sys.argv[4]   # e.g. chr12, chrX, etc
LABEL_SOURCE = sys.argv[5]  # e.g. histonemarks, tfitlatest
RANDOM_LABELS = False
if len(sys.argv) > 6:
    RANDOM_LABELS = True
print("Learning rate: %f" % LR)
print("Hidden Layer Size: %d" % HIDDEN_SIZE)
DATA_DIR = '/scratch/Users/igtr4848/atac_peak_features'


###############################################################################
## Load Data
print("Loading data...")
data = pd.read_pickle('%s/combined_dataset_union_%s.pkl' % (DATA_DIR, LABEL_SOURCE))
#data = pd.read_pickle('%s/combined_dataset_%s.pkl' % (DATA_DIR, LABEL_SOURCE))

print("TESTING ON %s - %s" % (TEST_SAMPLE, TEST_CHROM))

train_set = None
test_set = None
conditions_filename = ''
conditions_label = ''
if TEST_CHROM == 'all': # TEST_SAMPLE can't be "all"
    if TEST_SAMPLE == 'all':
        print('Invalid scenario')
        exit(0)

    train_set = data[(data['sample'] != TEST_SAMPLE)]
    test_set = data[(data['sample'] == TEST_SAMPLE)]
    conditions_filename = TEST_SAMPLE
    conditions_label = TEST_SAMPLE
else:
    if TEST_SAMPLE == 'all':
        train_set = data[(data.chrom != TEST_CHROM)]
        test_set = data[(data.chrom == TEST_CHROM)]
        conditions_filename = TEST_CHROM
        conditions_label = TEST_CHROM
    else:
        train_set = data[(data['chrom'] != TEST_CHROM) & (data['sample'] != TEST_SAMPLE)]
        test_set = data[(data['chrom'] == TEST_CHROM) & (data['sample'] == TEST_SAMPLE)]
        conditions_filename = "%s-%s" % (TEST_CHROM, TEST_SAMPLE)
        conditions_label = "%s / %s" % (TEST_CHROM, TEST_SAMPLE)

del data


###############################################################################
## Build Features
sigs = FeatureUnion([
    ('sigs', Pipeline([
        ('feature', ItemSelector('signal_features')),
        ('reshape', Reshape()),
    ])),
])

# sequence
S = set(train_set['sequence'])
T = Tokenizer(char_level=True, oov_token='*')
T.fit_on_texts(S)
seqs = FeatureUnion([
    ('seqs', Pipeline([
        ('feature', ItemSelector('sequence')),
        ('to_index', Seq2Ind(T))
    ]))
])


sigs.fit(train_set)
seqs.fit(train_set)


###############################################################################
## Build Classifier Pipeline
# which hyperparameters?
# rnn_hidden_size, learning_rate (let embedding_dim=108)

callbacks_list = [ReduceLROnPlateau(patience=3),
                  EarlyStopping(patience=4, restore_best_weights=True)]

pipeline = Pipeline([
    ('clf', KerasClassifier(build_fn=RNNHybridClassifier,
                            vocab_size=len(T.word_index),
                            input_length=1000,
                            embedding_dim=108,
                            rnn_hidden_size=HIDDEN_SIZE,
                            lr=LR,
                            batch_size=128,
                            epochs=40,
                            validation_split=0.1,
                            callbacks=callbacks_list,
                            nr_gpus=NR_GPUS,
                            verbose=1,))
])


low = 1*(train_set.mean_nr_nascent_reads < 0.001).values
high = 1*(train_set.mean_nr_nascent_reads > 0.01).values
usable = low + high
usable_idx = np.where(usable==1)[0]
train_set = train_set.iloc[usable_idx]
y_labels = 1*(train_set.mean_nr_nascent_reads > 0.01).values

my_y_test = 1*(test_set.mean_nr_nascent_reads > 0.01).values


if RANDOM_LABELS:
    #y_labels = np.random.randint(low=2, size=len(train_set))
    y_labels = np.random.choice(2, len(train_set), p=[0.36, 0.64])
###############################################################################
## Fit classifier and evaluate
pipeline.fit(X=[sigs.transform(train_set), seqs.transform(train_set)],
             y=y_labels)

## Evaluate performance on the test set using the best model
y_pred = pipeline.predict([sigs.transform(test_set), seqs.transform(test_set)])

if RANDOM_LABELS:
    my_y_test = np.random.choice(2, len(test_set), p=[0.36, 0.64])
coordinates = test_set[['chrom','start','end']].values
target_names = ['low txn', 'high txn']

print("-------------------------------------------------------------------")
if RANDOM_LABELS:
    print("RNN (Hybrid model) Confusion Matrix (RANDOM LABELS):")
else:
    print("RNN (Hybrid model) Confusion Matrix:")
print(confusion_matrix(my_y_test, y_pred))
print("-------------------------------------------------------------------")
print("F1-score (macro): %.3f" % f1_score(my_y_test, y_pred, average='macro'))
print("F1-score (micro): %.3f" % f1_score(my_y_test, y_pred, average='micro'))
print("F1-score (weighted): %.3f\n" % f1_score(my_y_test, y_pred, average='weighted'))
print("classification_report:")
print(classification_report(my_y_test, y_pred, target_names=target_names))
print("Accuracy: {}".format(accuracy_score(my_y_test, y_pred)))

## Plot ROC and report AUC
probas = pipeline.predict_proba([sigs.transform(test_set), seqs.transform(test_set)])
fpr, tpr, thresholds = roc_curve(my_y_test, probas[:, 1])
print("AUC: {}".format(auc(fpr, tpr)))
plt.plot(fpr, tpr);
if RANDOM_LABELS:
    plt.title('Receiver Operating Characteristic \n testing on %s (random labels,  AUC = %.3f)' % (conditions_label, auc(fpr, tpr)))
else:
    plt.title('Receiver Operating Characteristic \n testing on %s  (AUC = %.3f)' % (conditions_label, auc(fpr, tpr)))
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.axis('equal')
plt.xlim((0, 1));
plt.ylim((0, 1));
plt.plot(range(2), range(2), '--');
if RANDOM_LABELS:
    plt.savefig("%s/%s-FPFR-curve-random.png" % (DATA_DIR, conditions_filename), dpi=600)
else:
    plt.savefig("%s/%s-FPFR-curve.png" % (DATA_DIR, conditions_filename), dpi=600)

print('===========================================')
print('RF AUC: %.3f' % auc(fpr, tpr))
print('RF F1: %.3f' % f1_score(my_y_test, y_pred, average='weighted'))
print('RF Accuracy: %.3f' % accuracy_score(my_y_test, y_pred))
print('===========================================')


if not RANDOM_LABELS:
    model = pipeline.named_steps['clf'].model
    model_json = model.to_json()
    with open("%s/%s-hybrid-RNN-model.json" % (DATA_DIR, conditions_filename), "w") as json_file:
        json_file.write(model_json)
        # serialize weights to HDF5
        model.save_weights("%s/%s-hybrid-RNN-model-weights.h5" % (DATA_DIR, conditions_filename))
        print("Saved %s model to disk" % conditions_filename)

    false_positives = []
    false_negatives = []
    true_positives = []
    true_negatives = []
    for i in range(len(y_pred)):
        if my_y_test[i] == 0 and y_pred[i] == 1:
            false_positives.append(coordinates[i])
        elif my_y_test[i] == 1 and y_pred[i] == 0:
            false_negatives.append(coordinates[i])
        elif my_y_test[i] == 0:
            true_negatives.append(coordinates[i])
        else:
            true_positives.append(coordinates[i])

    with open("%s/%s-vs-%s_RNN-hybrid_false_positives.bed" % (DATA_DIR, conditions_filename, LABEL_SOURCE), 'w') as fp_file:
        for peak in false_positives:
            fp_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))
    with open("%s/%s-vs-%s_RNN-hybrid_false_negatives.bed" % (DATA_DIR, conditions_filename, LABEL_SOURCE), 'w') as fn_file:
        for peak in false_negatives:
            fn_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))
    with open("%s/%s-vs-%s_RNN-hybrid_true_positives.bed" % (DATA_DIR, conditions_filename, LABEL_SOURCE), 'w') as tp_file:
        for peak in true_positives:
            tp_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))
    with open("%s/%s-vs-%s_RNN-hybrid_true_negatives.bed" % (DATA_DIR, conditions_filename, LABEL_SOURCE), 'w') as tn_file:
        for peak in true_negatives:
            tn_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))

