import pandas as pd
import numpy as np
import multiprocessing
import matplotlib as mpl
import sys

# to prevent display weirdness when running in Pando:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


NR_THREADS = 24
EVALUATION_RADIUS = 500 # in bp
WHAT_TO_PLOT = sys.argv[1]

DATA_DIR = '/scratch/Users/igtr4848/atac_peak_features'
BEDGRAPH_DIR = '/scratch/Users/igtr4848/test/ATAC_peak_clf/mapped/rcc_bedgraphs'
PREFIX = 'rawnascent'


data = pd.read_pickle('%s/combined_dataset_union_%s.pkl' % (DATA_DIR, PREFIX))
data = data[data['sample'] == 'k562']

sample = 'k562'

def get_metapeaks(chromosome):
    tp_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    tp_neg_signal = np.zeros(2*EVALUATION_RADIUS)
    fp_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    fp_neg_signal = np.zeros(2*EVALUATION_RADIUS)
    tn_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    tn_neg_signal = np.zeros(2*EVALUATION_RADIUS)
    fn_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    fn_neg_signal = np.zeros(2*EVALUATION_RADIUS)

    pos_df = pd.read_csv('%s/%s.%s.pos.bedGraph' % (BEDGRAPH_DIR, sample, chromosome), header=None,
                          sep="\t", na_filter=False, usecols=[0, 1, 2, 3],
                          names=['chrom', 'start', 'end', 'reads'], \
                          dtype={'chrom':'str', 'start':'int', 'end':'int', 'reads':'float'})

    neg_df = pd.read_csv('%s/%s.%s.neg.bedGraph' % (BEDGRAPH_DIR, sample, chromosome), header=None,
                          sep="\t", na_filter=False, usecols=[0, 1, 2, 3],
                          names=['chrom', 'start', 'end', 'reads'], \
                          dtype={'chrom':'str', 'start':'int', 'end':'int', 'reads':'float'})

    if WHAT_TO_PLOT == 'tp':
        tp_file = "%s/%s-vs-rawnascent_RNN-hybrid_true_positives.bed" % (DATA_DIR, sample)

        tp_df = pd.read_csv(tp_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                            names=['chrom', 'start', 'end'],
                            dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
        tp_df = tp_df[tp_df['chrom'] == chromosome]

        for tp_region in tp_df.itertuples():
            # Get the number of reads for every nucleotide within the peak boundaries (both strands)
            tp_region_midpoint = tp_region.start + (tp_region.end - tp_region.start)/2
            window_start = tp_region_midpoint - EVALUATION_RADIUS
            window_end = tp_region_midpoint + EVALUATION_RADIUS
            pos_window_reads = np.zeros(2*EVALUATION_RADIUS)
            pos_signal_iter = pos_df.query('chrom == "%s" and start >= %d and end <= %d' % (tp_region.chrom, window_start, window_end)).itertuples()
            for coverage in pos_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    pos_window_reads[position] = coverage.reads
            tp_pos_signal = tp_pos_signal + pos_window_reads

            neg_window_reads = np.zeros(2*EVALUATION_RADIUS)
            neg_signal_iter = neg_df.query('chrom == "%s" and start >= %d and end <= %d' % (tp_region.chrom, window_start, window_end)).itertuples()
            for coverage in neg_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    neg_window_reads[position] = coverage.reads
            tp_neg_signal = tp_neg_signal + neg_window_reads

    elif WHAT_TO_PLOT == 'fp':
        fp_file = "%s/%s-vs-rawnascent_RNN-hybrid_false_positives.bed" % (DATA_DIR, sample)

        fp_df = pd.read_csv(fp_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                            names=['chrom', 'start', 'end'],
                            dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
        fp_df = fp_df[fp_df['chrom'] == chromosome]

        for fp_region in fp_df.itertuples():
            # Get the number of reads for every nucleotide within the peak boundaries (both strands)
            fp_region_midpoint = fp_region.start + (fp_region.end - fp_region.start)/2
            window_start = fp_region_midpoint - EVALUATION_RADIUS
            window_end = fp_region_midpoint + EVALUATION_RADIUS
            pos_window_reads = np.zeros(2*EVALUATION_RADIUS)
            pos_signal_iter = pos_df.query('chrom == "%s" and start >= %d and end <= %d' % (fp_region.chrom, window_start, window_end)).itertuples()
            for coverage in pos_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    pos_window_reads[position] = coverage.reads
            fp_pos_signal = fp_pos_signal + pos_window_reads

            neg_window_reads = np.zeros(2*EVALUATION_RADIUS)
            neg_signal_iter = neg_df.query('chrom == "%s" and start >= %d and end <= %d' % (fp_region.chrom, window_start, window_end)).itertuples()
            for coverage in neg_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    neg_window_reads[position] = coverage.reads
            fp_neg_signal = fp_neg_signal + neg_window_reads


    elif WHAT_TO_PLOT == 'tn':
        tn_file = "%s/%s-vs-rawnascent_RNN-hybrid_true_negatives.bed" % (DATA_DIR, sample)

        tn_df = pd.read_csv(tn_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                            names=['chrom', 'start', 'end'],
                            dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
        tn_df = tn_df[tn_df['chrom'] == chromosome]

        for tn_region in tn_df.itertuples():
            # Get the number of reads for every nucleotide within the peak boundaries (both strands)
            tn_region_midpoint = tn_region.start + (tn_region.end - tn_region.start)/2
            window_start = tn_region_midpoint - EVALUATION_RADIUS
            window_end = tn_region_midpoint + EVALUATION_RADIUS
            pos_window_reads = np.zeros(2*EVALUATION_RADIUS)
            pos_signal_iter = pos_df.query('chrom == "%s" and start >= %d and end <= %d' % (tn_region.chrom, window_start, window_end)).itertuples()
            for coverage in pos_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    pos_window_reads[position] = coverage.reads
            tn_pos_signal = tn_pos_signal + pos_window_reads

            neg_window_reads = np.zeros(2*EVALUATION_RADIUS)
            neg_signal_iter = neg_df.query('chrom == "%s" and start >= %d and end <= %d' % (tn_region.chrom, window_start, window_end)).itertuples()
            for coverage in neg_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    neg_window_reads[position] = coverage.reads
            tn_neg_signal = tn_neg_signal + neg_window_reads


    elif WHAT_TO_PLOT == 'fn':
        fn_file = "%s/%s-vs-rawnascent_RNN-hybrid_false_negatives.bed" % (DATA_DIR, sample)

        fn_df = pd.read_csv(fn_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                            names=['chrom', 'start', 'end'],
                            dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
        fn_df = fn_df[fn_df['chrom'] == chromosome]

        for fn_region in fn_df.itertuples():
            # Get the number of reads for every nucleotide within the peak boundaries (both strands)
            fn_region_midpoint = fn_region.start + (fn_region.end - fn_region.start)/2
            window_start = fn_region_midpoint - EVALUATION_RADIUS
            window_end = fn_region_midpoint + EVALUATION_RADIUS
            pos_window_reads = np.zeros(2*EVALUATION_RADIUS)
            pos_signal_iter = pos_df.query('chrom == "%s" and start >= %d and end <= %d' % (fn_region.chrom, window_start, window_end)).itertuples()
            for coverage in pos_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    pos_window_reads[position] = coverage.reads
            fn_pos_signal = fn_pos_signal + pos_window_reads

            neg_window_reads = np.zeros(2*EVALUATION_RADIUS)
            neg_signal_iter = neg_df.query('chrom == "%s" and start >= %d and end <= %d' % (fn_region.chrom, window_start, window_end)).itertuples()
            for coverage in neg_signal_iter:
                feat_start = coverage.start - window_start + 1
                feat_end = coverage.end - window_start
                for position in range(int(feat_start), int(feat_end)):
                    neg_window_reads[position] = coverage.reads
            fn_neg_signal = fn_neg_signal + neg_window_reads


    print("Done with chromosome %s" % chromosome)
    # return the partial sum of reads across this sample
    return tp_pos_signal, tp_neg_signal, fp_pos_signal, fp_neg_signal, tn_pos_signal, tn_neg_signal, fn_pos_signal, fn_neg_signal
          



if __name__=='__main__':

    CHROMOSOMES = [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', \
                    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', \
                    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', \
                    'chrX', 'chrY']


    pool = multiprocessing.Pool(NR_THREADS)
    final_result = pool.map(get_metapeaks, CHROMOSOMES)
    pool.close()
    pool.join()
    final_result = np.array(final_result)
    aggregated = final_result.sum(axis=0)

    my_dpi = 200

    if WHAT_TO_PLOT == 'tp':
        np.save('%s/k562-multi_nascent_TP_pos.npy' % DATA_DIR, aggregated[0])
        np.save('%s/k562-multi_nascent_TP_neg.npy' % DATA_DIR, aggregated[1])
        plt.clf()
        plt.rcParams.update({'font.size': 30})
        fig = plt.figure(figsize=(3000/my_dpi, 2500/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[0], color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[1], color='red')
        plt.title('Nascent transcription signal meta-peak \n from K562 true positives')
        plt.ylim((-850,850))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total # normalized reads')
        plt.tight_layout()
        plt.savefig("%s/k562-multi_TP_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'fp':
        np.save('%s/k562-multi_nascent_FP_pos.npy' % DATA_DIR, aggregated[2])
        np.save('%s/k562-multi_nascent_FP_neg.npy' % DATA_DIR, aggregated[3])
        plt.clf()
        plt.rcParams.update({'font.size': 30})
        fig = plt.figure(figsize=(3000/my_dpi, 2500/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[2], color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[3], color='red')
        plt.title('Nascent transcription signal meta-peak \n from K562 false positives')
        plt.ylim((-850,850))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total # normalized reads')
        plt.tight_layout()
        plt.savefig("%s/k562-multi_FP_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'tn':
        np.save('%s/k562-multi_nascent_TN_pos.npy' % DATA_DIR, aggregated[4])
        np.save('%s/k562-multi_nascent_TN_neg.npy' % DATA_DIR, aggregated[5])
        plt.clf()
        plt.rcParams.update({'font.size': 30})
        fig = plt.figure(figsize=(3000/my_dpi, 2500/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[4], color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[5], color='red')
        plt.title('Nascent transcription signal meta-peak \n from K562 true negatives')
        plt.ylim((-850,850))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total # normalized reads')
        plt.tight_layout()
        plt.savefig("%s/k562-multi_TN_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'fn':
        np.save('%s/k562-multi_nascent_FN_pos.npy' % DATA_DIR, aggregated[6])
        np.save('%s/k562-multi_nascent_FN_neg.npy' % DATA_DIR, aggregated[7])
        plt.clf()
        plt.rcParams.update({'font.size': 30})
        fig = plt.figure(figsize=(3000/my_dpi, 2500/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[6], color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), aggregated[7], color='red')
        plt.title('Nascent transcription signal meta-peak \n from K562 false negatives')
        plt.ylim((-850,850))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total # normalized reads')
        plt.tight_layout()
        plt.savefig("%s/k562-multi_FN_nascent_metapeak.png" % DATA_DIR)

