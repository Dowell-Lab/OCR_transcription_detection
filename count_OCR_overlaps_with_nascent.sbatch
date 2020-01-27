#!/bin/sh

cd ${DATADIR}

bedtools intersect -a ATAC/macs2/${1}_qe-6_dedup_peaks_clean.narrowPeak \
                   -b Nascent/fstitch_annotations/${1}_merged_fstitch_short_bidirs.bed \
                      Nascent/nf_tfit_annotations/${1}_merged_tfit_bidirs.bed \
                   -c >  ATAC/macs2/${1}_qe-6_dedup_peaks_clean_with_counts.narrowPeak
