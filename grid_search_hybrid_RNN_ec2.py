import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from tensorflow.keras.callbacks import ReduceLROnPlateau, EarlyStopping
#from keras import models
from sklearn.pipeline import FeatureUnion, Pipeline
from FeatureExtraction import ItemSelector, Reshape, Seq2Ind
from rnn_classifier import RNNHybridClassifier
from sklearn.metrics import confusion_matrix, f1_score, classification_report, roc_curve, auc, accuracy_score
from time import time
from keras.optimizers import Adam
import logging


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

log_filename = '/home/ec2-user/files/hyperparam_search_log.txt'
logging.basicConfig(filename=log_filename, level=logging.DEBUG)

train_samples = [
    'A549',
    'GM12878',
    'H1', 
    'HeLa',
    'K562',
    'LNCaP',
    'MCF7',
    'THP1'
    ]

###############################################################################
## data splits
val_samples = ['HCT116']

chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY']

split_index = 11 # ie. chr12
train_chroms = chromosomes[:split_index]
val_chroms = chromosomes[split_index:]

###############################################################################
## Load/Split Data
logging.debug("\n\nLoading data...")
data = pd.read_pickle('./data/combined_dataset_union_fstitchtfit.pkl')

train_set = data[data['sample'].isin(train_samples) & 
                 data['chrom'].isin(train_chroms)]

val_set = data[data['sample'].isin(val_samples) & 
               data['chrom'].isin(val_chroms)]

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

X_train_sigs = sigs.transform(train_set)
X_train_seqs = seqs.transform(train_set)
X_val_sigs = sigs.transform(val_set)
X_val_seqs = seqs.transform(val_set)
y_train = train_set['ovlp_txn'].values
y_val = val_set['ovlp_txn'].values

val_accs = []
params = []

for embedding_dim in [15, 50, 100]:
    for dropout in [0.1, 0.2, 0.3]:
        for hidden_size in [100, 200, 350, 500]:
            for learning_rate in [0.001, 0.0005, 0.0001]:


                ###############################################################################
                ## Build Classifier Pipeline

                callbacks_list = [ReduceLROnPlateau(patience=3),
                                  EarlyStopping(patience=4, restore_best_weights=True)]

                clf = KerasClassifier(build_fn=RNNHybridClassifier,
                                      vocab_size=len(T.word_index),
                                      input_length=1000,
                                      embedding_dim=embedding_dim,
                                      rnn_hidden_size=hidden_size,
                                      lr=learning_rate,
                                      dropout=dropout,
                                      mom=0.8932)

                ###############################################################################
                ## Fit classifier and evaluate
                start_time = time()
                hist = clf.fit(x=[X_train_sigs, X_train_seqs], y=y_train,
                               validation_data=([X_val_sigs, X_val_seqs], y_val),
                               batch_size=128,
                               callbacks=callbacks_list,
                               verbose=1,
                               shuffle=True,
                               epochs=25,
                              )

                elapsed_time = time() - start_time
                logging.debug("\n\nTotal training time: %.2f min" % (elapsed_time / 60))

                val_accuracy = hist.history['val_acc'][-1]
                logging.debug("-----------------------------------------------------------------------------")
                logging.debug("embedding_dim=%d , dropout=%1f , hidden_size=%d , learning_rate=%.4f" % (embedding_dim, dropout, hidden_size, learning_rate))
                logging.debug("validation accuracy = %.4f" % val_accuracy)
                logging.debug("-----------------------------------------------------------------------------")

                val_accs.append(val_accuracy)
                params.append("embedding_dim=%d , dropout=%1f , hidden_size=%d , learning_rate=%.4f" % (embedding_dim, dropout, hidden_size, learning_rate))


logging.debug("\n\nBest performance: %.4f" % max(val_accs))
logging.debug("Optimal hyperparams: %s" % params[np.argmax(val_accs)])
