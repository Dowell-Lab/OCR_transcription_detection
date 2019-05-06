#!/bin/sh

cd /scratch/Shares/dowell/processed_atacseq/macs2

# intersect all bed files
bedtools intersect -wa -a gm12878_peaks_union.broadPeak_ovlp-nascent.bed -b h1_peaks_union.broadPeak_ovlp-nascent.bed > tmp1.bed
bedtools intersect -wa -a tmp1.bed -b hct116_peaks_union.broadPeak_ovlp-nascent.bed > tmp2.bed
bedtools intersect -wa -a tmp2.bed -b k562_peaks_union.broadPeak_ovlp-nascent.bed > tmp3.bed
bedtools intersect -wa -a tmp3.bed -b cd4pos_peaks_union.broadPeak_ovlp-nascent.bed > tmp4.bed
bedtools intersect -wa -a tmp4.bed -b jurkat_peaks_union.broadPeak_ovlp-nascent.bed > tmp5.bed
bedtools intersect -wa -a tmp5.bed -b lncap_peaks_union.broadPeak_ovlp-nascent.bed > common.dupl.bed

# remove duplicate peak entries (looks like there's a bug in bedtools)
awk '!a[$0]++' common.dupl.bed > common_peaks_to_all_cell_types.bed

wc -l common_peaks_to_all_cell_types.bed

mv common_peaks_to_all_cell_types.bed /scratch/Users/igtr4848/atac_peak_features
