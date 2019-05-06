#!/bin/bash

cd /scratch/Shares/dowell/processed_atacseq/macs2

bedtools intersect -wa -a gm12878_peaks_union.broadPeak_ovlp-nascent.bed -b h1_peaks_union.broadPeak_ovlp-nascent.bed > foo
bedtools merge -i foo > int1.bed
wc -l int1.bed

bedtools intersect -wa -a int1.bed -b hct116_peaks_union.broadPeak_ovlp-nascent.bed > foo
bedtools merge -i foo > int2.bed
wc -l int2.bed

bedtools intersect -wa -a int2.bed -b k562_peaks_union.broadPeak_ovlp-nascent.bed > foo
bedtools merge -i foo > int3.bed
wc -l int3.bed

bedtools intersect -wa -a int3.bed -b cd4pos_peaks_union.broadPeak_ovlp-nascent.bed > foo
bedtools merge -i foo > int4.bed
wc -l int4.bed

bedtools intersect -wa -a int4.bed -b jurkat_peaks_union.broadPeak_ovlp-nascent.bed > foo
bedtools merge -i foo > int5.bed
wc -l int5.bed

bedtools intersect -wa -a int5.bed -b lncap_peaks_union.broadPeak_ovlp-nascent.bed > foo
bedtools merge -i foo > peaks_intersecting_all_samples.bed


echo "Total number of peaks intersecting all samples:"
wc -l peaks_intersecting_all_samples.bed

bedtools intersect -v -a gm12878_peaks_union.broadPeak_ovlp-nascent.bed -b peaks_intersecting_all_samples.bed > gm12878_peaks_union_no-common-peaks.broadPeak_ovlp-nascent.bed
bedtools intersect -v -a h1_peaks_union.broadPeak_ovlp-nascent.bed -b peaks_intersecting_all_samples.bed > h1_peaks_union_no-common-peaks.broadPeak_ovlp-nascent.bed
bedtools intersect -v -a hct116_peaks_union.broadPeak_ovlp-nascent.bed -b peaks_intersecting_all_samples.bed > hct116_peaks_union_no-common-peaks.broadPeak_ovlp-nascent.bed
bedtools intersect -v -a k562_peaks_union.broadPeak_ovlp-nascent.bed -b peaks_intersecting_all_samples.bed > k562_peaks_union_no-common-peaks.broadPeak_ovlp-nascent.bed
bedtools intersect -v -a cd4pos_peaks_union.broadPeak_ovlp-nascent.bed -b peaks_intersecting_all_samples.bed > cd4pos_peaks_union_no-common-peaks.broadPeak_ovlp-nascent.bed
bedtools intersect -v -a jurkat_peaks_union.broadPeak_ovlp-nascent.bed -b peaks_intersecting_all_samples.bed > jurkat_peaks_union_no-common-peaks.broadPeak_ovlp-nascent.bed
bedtools intersect -v -a lncap_peaks_union.broadPeak_ovlp-nascent.bed -b peaks_intersecting_all_samples.bed > lncap_peaks_union_no-common-peaks.broadPeak_ovlp-nascent.bed

#echo "Common peaks intersecting SRR1822165 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822165_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822165 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822165_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822165 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822165-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822165 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822165-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822165 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822165-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822165 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822165-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR1822166 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822166_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822166 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822166_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822166 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822166-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822166 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822166-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822166 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822166-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822166 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822166-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR1822167 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822167_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822167 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822167_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822167 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822167-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822167 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822167-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822167 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822167-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822167 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822167-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR1822168 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822168_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822168 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822168_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822168 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822168-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822168 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822168-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822168 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822168-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR1822168 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR1822168-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR5007258 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007258_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007258 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007258_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007258 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007258-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007258 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007258-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007258 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007258-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007258 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007258-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR5007259 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007259_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007259 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007259_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007259 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007259-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007259 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007259-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007259 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007259-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5007259 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5007259-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR5876158 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876158_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876158 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876158_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876158 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876158-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876158 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876158-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876158 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876158-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876158 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876158-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR5876159 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876159_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876159 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876159_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876159 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876159-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876159 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876159-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876159 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876159-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5876159 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5876159-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
#echo "Common peaks intersecting SRR5128074 hybrid RNN TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5128074_RNN-hybrid_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5128074 hybrid RNN FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5128074_RNN-hybrid_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5128074 attributes SVM TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5128074-vs-histonemarksandtfit-attr_SVM_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5128074 attributes SVM FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5128074-vs-histonemarksandtfit-attr_SVM_false_positives.bed | wc -l
#echo "Common peaks intersecting SRR5128074 signal RF TP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5128074-vs-histonemarksandtfit-sig_RF_true_positives.bed | wc -l
#echo "Common peaks intersecting SRR5128074 signal RF FP: "
#bedtools intersect -wa -a peaks_intersecting_all_samples.bed -b SRR5128074-vs-histonemarksandtfit-sig_RF_false_positives.bed | wc -l
#
