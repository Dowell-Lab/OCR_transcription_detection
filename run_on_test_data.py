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
TEST_SAMPLE = sys.argv[3]   # Sample prefix of the processed ATAC-seq sample to test (e.g. AllenGROnucnutlin)
LABEL_SOURCE = 'rawnascent'
print("Learning rate: %f" % LR)
print("Hidden Layer Size: %d" % HIDDEN_SIZE)
DATA_DIR = '/scratch/Users/igtr4848/atac_peak_features'

SUFFIX = '0.001-to-0.01'

###############################################################################
## Load Data
print("Loading data...")
train_set = pd.read_pickle('%s/combined_dataset_union_%s.pkl' % (DATA_DIR, LABEL_SOURCE))
test_set = pd.read_pickle('%s/%s-vs-rawnascent_from-raw-nascent-reads.pk' % (DATA_DIR, TEST_SAMPLE))

print("TESTING ON %s" % TEST_SAMPLE)


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


low = 1*(train_set.mean_nr_nascent_reads < 0.0005).values
high = 1*(train_set.mean_nr_nascent_reads > 0.03).values
usable = low + high
usable_idx = np.where(usable==1)[0]
train_set = train_set.iloc[usable_idx]

low = 1*(test_set.mean_nr_nascent_reads < 0.0005).values
high = 1*(test_set.mean_nr_nascent_reads > 0.03).values
usable = low + high
usable_idx = np.where(usable==1)[0]
test_set = test_set.iloc[usable_idx]

y_labels = 1*(train_set.mean_nr_nascent_reads > 0.03).values
my_y_test = 1*(test_set.mean_nr_nascent_reads > 0.03).values

###############################################################################
## Fit classifier and evaluate
pipeline.fit(X=[sigs.transform(train_set), seqs.transform(train_set)],
             y=y_labels)

## Evaluate performance on the test set using the best model
y_pred = pipeline.predict([sigs.transform(test_set), seqs.transform(test_set)])

coordinates = test_set[['chrom','start','end']].values
target_names = ['non-functional binding', 'functional binding']

positives = []
for i in range(len(y_pred)):
    if y_pred[i] == 1:
        positives.append(coordinates[i])

with open("%s/%s-vs-%s_RNN-hybrid_predicted_positives_%s.bed" % (DATA_DIR, TEST_SAMPLE, LABEL_SOURCE, SUFFIX), 'w') as fp_file:
    for peak in positives:
        fp_file.write("%s\t%s\t%s\n" % (peak[0], peak[1], peak[2]))

print("-------------------------------------------------------------------")
print("RNN (Hybrid model) Confusion Matrix_from-hct116-only:")
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
plt.title('Receiver Operating Characteristic \n testing on %s  (AUC = %.3f)' % (TEST_SAMPLE, auc(fpr, tpr)))
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.axis('equal')
plt.xlim((0, 1));
plt.ylim((0, 1));
plt.plot(range(2), range(2), '--');
plt.savefig("%s/%s-FPFR-curve_%s.png" % (DATA_DIR, TEST_SAMPLE, SUFFIX), dpi=600)


model = pipeline.named_steps['clf'].model
model_json = model.to_json()
with open("%s/%s-hybrid-RNN-model_%s.json" % (DATA_DIR, TEST_SAMPLE, SUFFIX), "w") as json_file:
    json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights("%s/%s-hybrid-RNN-model-weights_%s.h5" % (DATA_DIR, TEST_SAMPLE, SUFFIX))
    print("Saved %s model to disk" % TEST_SAMPLE)

print("Done")

