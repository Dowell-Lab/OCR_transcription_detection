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

#EMBEDDING_DIM = 120
#DROPOUT = 0.1
#HIDDEN_SIZE = 134
#LEARNING_RATE = 0.1597
EMBEDDING_DIM = 50
DROPOUT = 0.1
HIDDEN_SIZE = 100
LEARNING_RATE = 0.0001


log_filename = '/home/ec2-user/files/rnn_training_log.txt'
logging.basicConfig(filename=log_filename, level=logging.DEBUG)


###############################################################################
## Load/Split Data
logging.debug("\n\nLoading data...")
data = pd.read_pickle('./data/combined_dataset_union_fstitchtfit.pkl')

val_set = data[data['sample'] == 'HCT116']

train_set = data[data['sample'] != 'HCT116']

logging.debug("%d train samples, validating with %d samples" % (len(train_set), len(val_set)))
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


###############################################################################
## Build Classifier Pipeline

callbacks_list = [ReduceLROnPlateau(patience=3),
                  EarlyStopping(patience=4, restore_best_weights=True)]

clf = KerasClassifier(build_fn=RNNHybridClassifier,
                      vocab_size=len(T.word_index),
                      input_length=1000,
                      embedding_dim=EMBEDDING_DIM,
                      rnn_hidden_size=HIDDEN_SIZE,
                      lr=LEARNING_RATE,
                      dropout=DROPOUT,
                      mom=0.8932)

###############################################################################
## Fit classifier and evaluate
start_time = time()
clf.fit(x=[X_train_sigs, X_train_seqs], y=y_train,
        validation_data=([X_val_sigs, X_val_seqs], y_val),
        batch_size=128,
        callbacks=callbacks_list,
        verbose=1,
        shuffle=True,
        epochs=25,
       )

tf.keras.models.save_model(clf.model, f"./data/saved_models/hybrid_all.h5")
#    saved_model = tf.keras.models.load_model(f"./data/saved_models/hybrid_{test_sample}.h5")
elapsed_time = time() - start_time
logging.debug("Total %s training time: %.2f min" % (test_sample, (elapsed_time / 60)))

