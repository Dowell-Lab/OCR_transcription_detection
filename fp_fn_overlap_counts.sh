#!/bin/bash

cd /scratch/Users/igtr4848/atac_peak_features

wc -l SRR5128074-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed
wc -l SRR5128074-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed
wc -l SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed
bedtools intersect -wa -b SRR5128074-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed -a SRR5128074-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed > rnn-attr.bed
wc -l rnn-attr.bed
bedtools intersect -wa -a rnn-attr.bed -b SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l
bedtools intersect -wa -b SRR5128074-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed -a SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l
bedtools intersect -wa -b SRR5128074-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed -a SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l



wc -l SRR1822165-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed
wc -l SRR1822165-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed
wc -l SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed
bedtools intersect -wa -b SRR1822165-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed -a SRR1822165-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed > rnn-attr.bed
wc -l rnn-attr.bed
bedtools intersect -wa -a rnn-attr.bed -b SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l
bedtools intersect -wa -b SRR1822165-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed -a SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l
bedtools intersect -wa -b SRR1822165-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed -a SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l



wc -l SRR5007258-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed
wc -l SRR5007258-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed
wc -l SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed
bedtools intersect -wa -b SRR5007258-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed -a SRR5007258-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed > rnn-attr.bed
wc -l rnn-attr.bed
bedtools intersect -wa -a rnn-attr.bed -b SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l
bedtools intersect -wa -b SRR5007258-vs-tfitfstitchhistmarks_RNN-hybrid_false_positives.bed -a SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l
bedtools intersect -wa -b SRR5007258-vs-tfitfstitchhistmarks-attr_SVM_false_positives.bed -a SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_positives.bed | wc -l



wc -l SRR5128074-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed
wc -l SRR5128074-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed
wc -l SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed
bedtools intersect -wa -b SRR5128074-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed -a SRR5128074-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed > rnn-attr.bed
wc -l rnn-attr.bed
bedtools intersect -wa -a rnn-attr.bed -b SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l
bedtools intersect -wa -b SRR5128074-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed -a SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l
bedtools intersect -wa -b SRR5128074-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed -a SRR5128074-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l



wc -l SRR1822165-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed
wc -l SRR1822165-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed
wc -l SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed
bedtools intersect -wa -b SRR1822165-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed -a SRR1822165-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed > rnn-attr.bed
wc -l rnn-attr.bed
bedtools intersect -wa -a rnn-attr.bed -b SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l
bedtools intersect -wa -b SRR1822165-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed -a SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l
bedtools intersect -wa -b SRR1822165-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed -a SRR1822165-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l



wc -l SRR5007258-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed
wc -l SRR5007258-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed
wc -l SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed
bedtools intersect -wa -b SRR5007258-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed -a SRR5007258-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed > rnn-attr.bed
wc -l rnn-attr.bed
bedtools intersect -wa -a rnn-attr.bed -b SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l
bedtools intersect -wa -b SRR5007258-vs-_tfitfstitchhistmarks_RNN-hybrid_false_negatives.bed -a SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l
bedtools intersect -wa -b SRR5007258-vs-tfitfstitchhistmarks-attr_SVM_false_negatives.bed -a SRR5007258-vs-tfitfstitchhistmarks-sig_RF_false_negatives.bed | wc -l



