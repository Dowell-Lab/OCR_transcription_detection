#!/bin/sh

cd ${ATACPEAKSDIR}

# intersect all bed files (except for the HCT116 cells, since those are only used for validation)
bedtools intersect -wa -a A549_qe-6_dedup_peaks_clean_with_counts.narrowPeak -b GM12878_qe-6_dedup_peaks_clean_with_counts.narrowPeak > tmp1.bed
bedtools intersect -wa -a tmp1.bed -b H1_qe-6_dedup_peaks_clean_with_counts.narrowPeak > tmp2.bed
bedtools intersect -wa -a tmp2.bed -b LNCaP_qe-6_dedup_peaks_clean_with_counts.narrowPeak > tmp3.bed
bedtools intersect -wa -a tmp3.bed -b HeLa_qe-6_dedup_peaks_clean_with_counts.narrowPeak > tmp4.bed
bedtools intersect -wa -a tmp4.bed -b K562_qe-6_dedup_peaks_clean_with_counts.narrowPeak > tmp5.bed
bedtools intersect -wa -a tmp5.bed -b MCF7_qe-6_dedup_peaks_clean_with_counts.narrowPeak > tmp6.bed
bedtools intersect -wa -a tmp6.bed -b THP1_qe-6_dedup_peaks_clean_with_counts.narrowPeak > common.dupl.bed

# remove duplicate peak entries (looks like there's a bug in bedtools)
awk '!a[$0]++' common.dupl.bed > common_peaks_to_all_cell_types.bed

wc -l common_peaks_to_all_cell_types.bed

