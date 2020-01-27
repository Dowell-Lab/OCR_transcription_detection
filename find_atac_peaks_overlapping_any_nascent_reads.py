import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import argparse
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import multiprocessing

import matplotlib as mpl 

mpl.use('Agg')
import matplotlib.pyplot as plt 
plt.ioff()

NR_THREADS = 24
SIGNAL_EVALUATION_RADIUS = 500
SEQUENCE_EVALUATION_RADIUS = 500
SEQTK_PATH = '/Users/igtr4848/seqtk/seqtk'

parser = argparse.ArgumentParser(description='This script gathers all features that can be used to classify ATAC-seq peaks')
parser.add_argument('-p', '--atac-peaks', dest='atac_peaks', help='Path to the ATAC-seq peaks file')
parser.add_argument('-n', '--nascent-bedgraph', dest='nascent_bedgraph', help='Path to the nascent bedgraph for a matching cell type and conditions')
parser.add_argument('-m', '--min-avg-reads', dest='min_read_threshold', help='Minimum average number of nascent transcription reads within the evaluation radius', required=False, default=0.01, type=float)
args = parser.parse_args()  # /scratch/Shares/dowell/nascent/hg38/hg38_refseq.bed


atac_peak_df = pd.read_csv(args.atac_peaks, header=None, \
                           sep="\t", na_filter=False, \
                           usecols=[0, 1, 2], \
                           names=['chrom', 'start', 'end'], \
                           dtype={'chrom':'str', 'start':'int', \
                                   'end':'int'} )

nascent_bed_df = pd.read_csv( args.nascent_bedgraph, header=None, \
                              sep="\t", na_filter=False, \
                              usecols=[0, 1, 2, 3], \
                              names=['chrom', 'start', 'end', 'reads'], \
                              dtype={'chrom':'str', 'start':'int', 'end':'int', 'reads':'float'} )

def get_nascent_coverage_for_peaks(current_chrom):
    print("Processing peaks from chromosome %s" % current_chrom)
    nascent_coverage = []
    peaks_with_nascent = []
    ovlp_iter = atac_peak_df[(atac_peak_df.chrom == current_chrom)].itertuples()

    for peak in ovlp_iter:
        peak_midpoint = peak.start + (peak.end - peak.start)/2
        region_start = peak_midpoint - SIGNAL_EVALUATION_RADIUS
        region_end = peak_midpoint + SIGNAL_EVALUATION_RADIUS
        peak_window_reads = np.zeros(2*SIGNAL_EVALUATION_RADIUS)

        # Get the number of reads for every nucleotide +/- the SIGNAL_EVALUATION_RADIUS
        wide_signal_iter = nascent_bed_df.query('chrom == "%s" and start >= %d and end <= %d' % (current_chrom, region_start, region_end)).itertuples()
        for coverage in wide_signal_iter:
            feat_start = coverage.start - region_start + 1
            feat_end = coverage.end - region_start
            for position in range(int(feat_start), int(feat_end)):
                peak_window_reads[position] = abs(coverage.reads)

        mean_nr_reads = np.mean(peak_window_reads)
        peaks_with_nascent.append([peak.chrom, peak.start, peak.end, mean_nr_reads])

    return(peaks_with_nascent)


if __name__=='__main__':

    CHROMOSOMES = [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', \
                    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', \
                    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', \
                    'chrX', 'chrY']

    pool = multiprocessing.Pool(NR_THREADS)
    results = pool.map(get_nascent_coverage_for_peaks, CHROMOSOMES)
    pool.close()
    pool.join()

    with open("%s_ovlp-nascent.bed" % args.atac_peaks, 'w') as filtered_peaks:
        for chrom in results:
            for peak in chrom:
                if len(peak) > 0:
                    filtered_peaks.write("%s\t%s\t%s\t%s\n" % (peak[0], peak[1], peak[2], peak[3]))

