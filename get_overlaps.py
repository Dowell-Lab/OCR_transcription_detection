import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import argparse
from argparse import ArgumentParser
from Bio import SeqIO
import numpy as np
import os
import pandas as pd
import pickle
import multiprocessing

NR_THREADS = 24
SIGNAL_EVALUATION_RADIUS = 500
SEQUENCE_EVALUATION_RADIUS = 500
SEQTK_PATH = '/Users/igtr4848/seqtk/seqtk'
SEQ_DUMP_DIR = '/scratch/Shares/dowell/itripodi/ATAC/peak_sequences'
OUTDIR = '/scratch/Shares/dowell/itripodi/ATAC/features'

parser = argparse.ArgumentParser(description='This script gathers all features that can be used to classify ATAC-seq peaks')
parser.add_argument('-x', '--prefix', dest='prefix', help='Custom prefix for all output files')
parser.add_argument('-p', '--atac-peaks', dest='atac_peaks', help='Path to the ATAC-seq peaks file')
parser.add_argument('-b', '--atac-bedgraph', dest='atac_bedgraph', help='Path to the ATAC-seq bedGraph file for the same sample of the given peaks file.')
parser.add_argument('-f', '--genome-fasta', dest='genome_fasta', help='Path to the FASTA file for the reference genome in use.')
parser.add_argument('-r', '--refseq-annotations', dest='refseq_annotations', help='Path to the RefSeq annotation file for the same reference genome these files were processed.')
parser.add_argument('-s', '--skip-seq-caching', dest='skip_sequence_caching', help='Skip the caching of peak sequences if this has already been processed.', required=False, action='store_true')
args = parser.parse_args()  # /scratch/Shares/dowell/nascent/hg38/hg38_refseq.bed


atac_peak_df = pd.read_csv(args.atac_peaks, header=None, \
                           sep="\t", na_filter=False, \
                           usecols=[0, 1, 2, 10], \
                           names=['chrom', 'start', 'end', 'ovlp_txn'], \
                           dtype={'chrom':'str', 'start':'int', \
                                  'end':'int', 'ovlp_txn':'int'} )

atac_bed_df = pd.read_csv(args.atac_bedgraph, header=None, \
                           sep="\t", na_filter=False, \
                           usecols=[0, 1, 2, 3], \
                           names=['chrom', 'start', 'end', 'reads'], \
                           dtype={'chrom':'str', 'start':'int', \
                                  'end':'int', 'reads':'float'} )

refseq_df = pd.read_csv(args.refseq_annotations, header=None, \
                           sep="\t", na_filter=False, \
                           usecols=[0, 1, 2, 5], \
                           names=['chrom', 'start', 'end', 'strand'], \
                           dtype={'chrom':'str', 'start':'int', \
                                  'end':'int', 'strand':'str'} )

def get_atac_features(current_chrom):
    print("Processing peaks from chromosome %s" % current_chrom)
    chrom_features = []
    last_peak_end = 0
    ovlp_iter = atac_peak_df[(atac_peak_df.chrom == current_chrom)].itertuples()

    for peak in ovlp_iter:
        signal_features = np.zeros(2 * SIGNAL_EVALUATION_RADIUS)
        peak_midpoint = peak.start + (peak.end - peak.start)/2
        region_start = peak_midpoint - SIGNAL_EVALUATION_RADIUS
        region_end = peak_midpoint + SIGNAL_EVALUATION_RADIUS
        peak_window_reads = np.zeros(peak.end - peak.start)

        # Get the number of reads for every nucleotide +/- the SIGNAL_EVALUATION_RADIUS
        #wide_signal_iter = atac_bed_df[(atac_bed_df['chrom'] == current_chrom) & (atac_bed_df['start'] >= region_start) & (atac_bed_df['end'] <= region_end)].itertuples()
        wide_signal_iter = atac_bed_df.query('chrom == "%s" and start >= %d and end <= %d' % (current_chrom, region_start, region_end)).itertuples()
        for coverage in wide_signal_iter:
            feat_start = coverage.start - region_start + 1
            feat_end = coverage.end - region_start
            for position in range(int(feat_start), int(feat_end)):
                signal_features[position] = coverage.reads

        # Get the mean number of reads for every nucleotide within the peak boundaries
        peak_signal_iter = atac_bed_df.query('chrom == "%s" and start >= %d and end <= %d' % (current_chrom, peak.start, peak.end)).itertuples()
        for coverage in peak_signal_iter:
            feat_start = coverage.start - peak.start + 1
            feat_end = coverage.end - peak.start
            for position in range(int(feat_start), int(feat_end)):
                peak_window_reads[position] = coverage.reads

        # Gather other peak attributes
        dist_from_last_peak = 0
        if last_peak_end > 0:
            dist_from_last_peak = peak.start - last_peak_end
            last_peak_end = peak.end
        peak_width = peak.end - peak.start

        # Check whether it overlaps a promoter
        promoter_overlaps = len(refseq_df.query('chrom == "%s" and ((strand == "+" and %d < start and %d > start) or (strand == "-" and %d < end and %d > end))' % (current_chrom, peak.start, peak.end, peak.start, peak.end)))

        # Get sequence-based features
        sequence_row = seq_df.query('chrom == "%s" and start < %d and end > %d' % (current_chrom, peak_midpoint, peak_midpoint))
        if len(sequence_row) > 0:
            sequence = sequence_row.iloc[0].seq
            gc_count = 0
            for nucleotide in sequence:
                if nucleotide == 'C' or nucleotide == 'G' or nucleotide == 'S':
                    gc_count += 1
            gc_ratio = float(gc_count) / (2 * SEQUENCE_EVALUATION_RADIUS)
        else:
            print("Could not find sequence information for peak %s:%s-%s ... skipping peak" % (current_chrom, peak.start, peak.end))
            continue

        overlaps_nascent_txn = 0
        if peak.ovlp_txn > 0:
            overlaps_nascent_txn = 1

        chrom_features.append(  { 'chrom': current_chrom, \
                                  'start': peak.start, \
                                  'end': peak.end, \
                                  'ovlp_txn': overlaps_nascent_txn, \
                                  'prom_ovlp': promoter_overlaps, \
                                  'width': peak_width, \
                                  'mean_nr_reads': np.mean(peak_window_reads), \
                                  'max_reads': np.max(peak_window_reads), \
                                  'min_reads': np.min(peak_window_reads), \
                                  'dist_from_last_peak': dist_from_last_peak, \
                                  'gc_ratio': gc_ratio,
                                  'sequence': sequence,
                                  'signal_features': signal_features } )
        last_peak_end = peak.end
    return np.array(chrom_features)


# easy way to create an empty pandas dataframe from
# https://stackoverflow.com/questions/36462257/create-empty-dataframe-in-pandas-specifying-column-types#36463086
def df_empty(columns, dtypes, index=None):
    assert len(columns)==len(dtypes)
    df = pd.DataFrame(index=index)
    for c,d in zip(columns, dtypes):
        df[c] = pd.Series(dtype=d)
    return df



if __name__=='__main__':

    CHROMOSOMES = [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', \
                    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', \
                    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', \
                    'chrX', 'chrY']

    atac_peaks_file = args.atac_peaks.split('/')[-1].split('.')[0]
    if args.skip_sequence_caching:  # don't waste cycles if this has already been cached
        seq_df = pd.read_pickle('%s/%s_seq_pandas_dataframe.pickle' % (SEQ_DUMP_DIR, atac_peaks_file))
    else:
        print('Gathering the actual sequence overlapping each peak:')
        # Generate a temp bed file with fixed windows for each peak
        tmp_peaks_filename = '%s/%s_tmp_fixed_peak_windows.bed' % (SEQ_DUMP_DIR, atac_peaks_file)
        with open(tmp_peaks_filename, 'w') as tmp_peaks:
            for peak in atac_peak_df.itertuples():
                if peak.chrom in CHROMOSOMES:
                    peak_midpoint = peak.start + (peak.end - peak.start)/2
                    tmp_peaks.write("%s\t%d\t%d\n" % (peak.chrom, \
                                                      peak_midpoint - SEQUENCE_EVALUATION_RADIUS, \
                                                      peak_midpoint + SEQUENCE_EVALUATION_RADIUS))
        # Gather the sequence for all those fixed-width regions
        sequence = os.popen('%s subseq %s %s > %s/%s.fa' % (SEQTK_PATH, args.genome_fasta, tmp_peaks_filename, SEQ_DUMP_DIR, atac_peaks_file)).read()
        fasta_sequences = SeqIO.parse(open('%s/%s.fa' % (SEQ_DUMP_DIR, atac_peaks_file)),'fasta')
        # Create a dataframe to access these more efficiently
        seq_df = df_empty(['chrom', 'start', 'end', 'seq'], [np.str, np.int, np.int, np.str])
        for seq in fasta_sequences:
            try:
                chrom, coordinates = seq.id.split(':')
                start = coordinates.split('-')[0]
                end = coordinates.split('-')[1]
                seq_df.loc[len(seq_df)] = [chrom, int(start), int(end), str(seq.seq)]
            except Exception as e:
                print("Something went sideways parsing the FASTA:")
                print(seq)
                print('-------------------------------------------')
                print(e)
                raise

        seq_df.to_pickle('%s/%s_seq_pandas_dataframe.pickle' % (SEQ_DUMP_DIR, atac_peaks_file))

    pool = multiprocessing.Pool(NR_THREADS)
    X_and_y_train = pool.map(get_atac_features, CHROMOSOMES)
    pool.close()
    pool.join()

    data = [j for sub in X_and_y_train for j in sub]
    data = pd.DataFrame(data)
    data['signal_features'] = data['signal_features'].map(lambda x: x.astype('float32'))
    data['ovlp_txn'] = data['ovlp_txn'].astype('int')
    data['mean_nr_reads'] = data['mean_nr_reads'].astype('float32')
    data['gc_ratio'] = data['gc_ratio'].astype('float32')
    data.to_pickle('%s/%s_from-fstitchtfit.pk' % (OUTDIR, args.prefix))

