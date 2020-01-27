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

EMBEDDING_DIM = 50
DROPOUT = 0.1
HIDDEN_SIZE = 100
LEARNING_RATE = 0.0001


log_filename = '/home/ec2-user/files/TSS-label_rnn_training_log.txt'
logging.basicConfig(filename=log_filename, level=logging.DEBUG)

###############################################################################
## data splits
val_samples = ['HCT116']

chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY']

split_index = 11 # ie. chr12
train_chroms = chromosomes[:split_index]
val_chroms = chromosomes[split_index:]


for test_index in range(8):
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
    test_sample = train_samples[test_index]
    logging.debug("\n\nTesting on %s" % test_sample)
    del train_samples[test_index]

    ###############################################################################
    ## Load/Split Data
    logging.debug("\n\nLoading data...")
    data = pd.read_pickle('./data/combined_dataset_union_fstitchtfit.pkl')

    train_set = data[data['sample'].isin(train_samples) & 
                     data['chrom'].isin(train_chroms)]

    val_set = data[data['sample'].isin(val_samples) & 
                   data['chrom'].isin(val_chroms)]

    test_set = data[(data['sample'] == test_sample) & 
                    data['chrom'].isin(val_chroms)]

    logging.debug("%d train samples, testing on %d samples" % (len(train_set), len(test_set)))
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
    X_test_sigs = sigs.transform(test_set)
    X_test_seqs = seqs.transform(test_set)
    y_train = train_set['prom_ovlp'].values
    y_val = val_set['ovlp_txn'].values
    y_test = test_set['ovlp_txn'].values


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

    tf.keras.models.save_model(clf.model, f"./data/saved_models/TSS-label_hybrid_{test_sample}.h5")
#    saved_model = tf.keras.models.load_model(f"./data/saved_models/TSS-label_hybrid_{test_sample}.h5")
    elapsed_time = time() - start_time
    logging.debug("Total %s training time: %.2f min" % (test_sample, (elapsed_time / 60)))

    y_prob = clf.model.predict(x=[X_test_sigs, X_test_seqs])
#    y_prob = saved_model.predict(x=[X_test_sigs, X_test_seqs])
    np.save('/home/ec2-user/files/TSS-label_%s_y_prob.npy' % test_sample, y_prob)
    #y_pred = y_prob.argmax(axis=-1)
    y_pred = [int(np.round(x)) for x in y_prob]
    np.save('/home/ec2-user/files/TSS-label_%s_y_test.npy' % test_sample, y_test)
    np.save('/home/ec2-user/files/TSS-label_%s_y_pred.npy' % test_sample, y_pred)
    elapsed_time = time() - start_time
    logging.debug("Total %s training and testing time: %.2f min" % (test_sample, (elapsed_time / 60)))
    target_names = ['No transcription', 'Transcription']

    logging.debug("-------------------------------------------------------------------")
    logging.debug("RNN (Hybrid model) %s Confusion Matrix:" % test_sample)
    logging.debug(confusion_matrix(y_test, y_pred))
    logging.debug("-------------------------------------------------------------------")
    logging.debug("F1-score (macro): %.3f" % f1_score(y_test, y_pred, average='macro'))
    logging.debug("F1-score (micro): %.3f" % f1_score(y_test, y_pred, average='micro'))
    logging.debug("F1-score (weighted): %.3f\n" % f1_score(y_test, y_pred, average='weighted'))
    logging.debug("classification_report:")
    logging.debug(classification_report(y_test, y_pred, target_names=target_names))
    logging.debug("Accuracy: {}".format(accuracy_score(y_test, y_pred)))

    ## Plot ROC and report AUC
    #probas = clf.predict_proba([sigs.transform(test_set), seqs.transform(test_set)])
    #np.save('/home/ec2-user/files/%s_probas.npy' % test_sample, probas)
    #fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
    fpr, tpr, thresholds = roc_curve(y_test, y_prob)
    np.save('/home/ec2-user/files/TSS-label_%s_fpr.npy' % test_sample, fpr)
    np.save('/home/ec2-user/files/TSS-label_%s_tpr.npy' % test_sample, tpr)
    np.save('/home/ec2-user/files/TSS-label_%s_thresholds.npy' % test_sample, thresholds)

    logging.debug("AUC: {}".format(auc(fpr, tpr)))
    np.save('/home/ec2-user/files/TSS-label_%s_auc.npy' % test_sample, auc)

    logging.debug('===========================================')
    logging.debug("Results for %s" % test_sample)
    logging.debug('RNN AUC: %.3f' % auc(fpr, tpr))
    logging.debug('RNN F1: %.3f' % f1_score(y_test, y_pred, average='weighted'))
    logging.debug('RNN Accuracy: %.3f' % accuracy_score(y_test, y_pred))
    logging.debug('===========================================')


    test_set_tss = test_set[test_set['prom_ovlp'] > 0]
    test_set_notss = test_set[test_set['prom_ovlp'] == 0]
    y_test_tss = test_set_tss['ovlp_txn'].values
    y_test_notss = test_set_notss['ovlp_txn'].values

    y_prob_tss = clf.model.predict(x=[sigs.transform(test_set_tss), seqs.transform(test_set_tss)])
    #y_prob_tss = saved_model.predict(x=[sigs.transform(test_set_tss), seqs.transform(test_set_tss)])
    #y_pred_tss = y_prob_tss.argmax(axis=-1)
    y_pred_tss = [int(np.round(x)) for x in y_prob_tss]

    logging.debug("-------------------------------------------------------------------")
    logging.debug("TSS RNN (Hybrid model) %s Confusion Matrix:" % test_sample)
    logging.debug(confusion_matrix(y_test_tss, y_pred_tss))
    logging.debug("-------------------------------------------------------------------")
    logging.debug("F1-score (weighted): %.3f\n" % f1_score(y_test_tss, y_pred_tss, average='weighted'))
    logging.debug("classification_report:")
    logging.debug(classification_report(y_test_tss, y_pred_tss, target_names=target_names))
    logging.debug("Accuracy: {}".format(accuracy_score(y_test_tss, y_pred_tss)))

    ## Plot ROC and report AUC
    #probas = clf.predict_proba([sigs.transform(test_set_tss), seqs.transform(test_set_tss)])
    #np.save('/home/ec2-user/files/tss_%s_probas.npy' % test_sample, probas)
    #fpr, tpr, thresholds = roc_curve(y_test_tss, probas[:, 1])
    fpr, tpr, thresholds = roc_curve(y_test_tss, y_prob_tss)
    np.save('/home/ec2-user/files/TSS-label_tss_%s_fpr.npy' % test_sample, fpr)
    np.save('/home/ec2-user/files/TSS-label_tss_%s_tpr.npy' % test_sample, tpr)
    np.save('/home/ec2-user/files/TSS-label_tss_%s_thresholds.npy' % test_sample, thresholds)
    np.save('/home/ec2-user/files/TSS-label_tss_%s_y_pred.npy' % test_sample, y_pred_tss)
    np.save('/home/ec2-user/files/TSS-label_tss_%s_y_prob.npy' % test_sample, y_prob_tss)
    np.save('/home/ec2-user/files/TSS-label_tss_%s_y_test.npy' % test_sample, y_test_tss)

    logging.debug("AUC: {}".format(auc(fpr, tpr)))
    np.save('/home/ec2-user/files/TSS-label_%s_auc.npy' % test_sample, auc(fpr, tpr))

    logging.debug('===========================================')
    logging.debug("Results for %s" % test_sample)
    logging.debug('TSS RNN AUC: %.3f' % auc(fpr, tpr))
    logging.debug('TSS RNN F1: %.3f' % f1_score(y_test_tss, y_pred_tss, average='weighted'))
    logging.debug('TSS RNN Accuracy: %.3f' % accuracy_score(y_test_tss, y_pred_tss))
    logging.debug('===========================================')


    y_prob_notss = clf.model.predict(x=[sigs.transform(test_set_notss), seqs.transform(test_set_notss)])
    #y_prob_notss = saved_model.predict(x=[sigs.transform(test_set_notss), seqs.transform(test_set_notss)])
    #y_pred_notss = y_prob_notss.argmax(axis=-1)
    y_pred_notss = [int(np.round(x)) for x in y_prob_notss]

    logging.debug("-------------------------------------------------------------------")
    logging.debug("Non-TSS RNN (Hybrid model) %s Confusion Matrix:" % test_sample)
    logging.debug(confusion_matrix(y_test_notss, y_pred_notss))
    logging.debug("-------------------------------------------------------------------")
    logging.debug("F1-score (weighted): %.3f\n" % f1_score(y_test_notss, y_pred_notss, average='weighted'))
    logging.debug("classification_report:")
    logging.debug(classification_report(y_test_notss, y_pred_notss, target_names=target_names))
    logging.debug("Accuracy: {}".format(accuracy_score(y_test_notss, y_pred_notss)))

    ## Plot ROC and report AUC
    #probas = clf.predict_proba([sigs.transform(test_set_notss), seqs.transform(test_set_notss)])
    #np.save('/home/ec2-user/files/notss_%s_probas.npy' % test_sample, probas)
    #fpr, tpr, thresholds = roc_curve(y_test_notss, probas[:, 1])
    fpr, tpr, thresholds = roc_curve(y_test_notss, y_pred_notss)
    np.save('/home/ec2-user/files/TSS-label_notss_%s_fpr.npy' % test_sample, fpr)
    np.save('/home/ec2-user/files/TSS-label_notss_%s_tpr.npy' % test_sample, tpr)
    np.save('/home/ec2-user/files/TSS-label_notss_%s_thresholds.npy' % test_sample, thresholds)
    np.save('/home/ec2-user/files/TSS-label_notss_%s_y_pred.npy' % test_sample, y_pred_notss)
    np.save('/home/ec2-user/files/TSS-label_notss_%s_y_prob.npy' % test_sample, y_prob_notss)
    np.save('/home/ec2-user/files/TSS-label_notss_%s_y_test.npy' % test_sample, y_test_notss)

    logging.debug("AUC: {}".format(auc(fpr, tpr)))
    np.save('/home/ec2-user/files/TSS-label_%s_auc.npy' % test_sample, auc(fpr, tpr))

    logging.debug('===========================================')
    logging.debug("Results for %s" % test_sample)
    logging.debug('Non-TSS RNN AUC: %.3f' % auc(fpr, tpr))
    logging.debug('Non-TSS RNN F1: %.3f' % f1_score(y_test_notss, y_pred_notss, average='weighted'))
    logging.debug('Non-TSS RNN Accuracy: %.3f' % accuracy_score(y_test_notss, y_pred_notss))
    logging.debug('===========================================')

