import pandas as pd
import numpy as np
import multiprocessing
import matplotlib as mpl
import sys

# to prevent display weirdness when running in Pando:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


NR_THREADS = 6
EVALUATION_RADIUS = 500 # in bp
WHAT_TO_PLOT = sys.argv[1]

DATA_DIR = './new_models/results/bed_files/'
BEDGRAPH_DIR = '/home/ignaciot/dowell_lab/rnn_paper_backup/nascent_bg'


def get_metapeaks(chromosome):
    tp_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    tp_neg_signal = np.zeros(2*EVALUATION_RADIUS)
    fp_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    fp_neg_signal = np.zeros(2*EVALUATION_RADIUS)
    tn_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    tn_neg_signal = np.zeros(2*EVALUATION_RADIUS)
    fn_pos_signal = np.zeros(2*EVALUATION_RADIUS)
    fn_neg_signal = np.zeros(2*EVALUATION_RADIUS)
    count = 0

    for sample in ['A549', 'GM12878', 'H1', 'HeLa', 'K562', 'LNCaP', 'MCF7', 'THP1']:
        pos_df = pd.read_csv('%s/%s.%s.pos.merged.clean.sorted.bedGraph' % (BEDGRAPH_DIR, sample, chromosome), header=None,
                            sep="\t", na_filter=False, usecols=[0, 1, 2, 3],
                            names=['chrom', 'start', 'end', 'reads'], \
                            dtype={'chrom':'str', 'start':'int', 'end':'int', 'reads':'float'})

        neg_df = pd.read_csv('%s/%s.%s.neg.merged.clean.sorted.bedGraph' % (BEDGRAPH_DIR, sample, chromosome), header=None,
                            sep="\t", na_filter=False, usecols=[0, 1, 2, 3],
                            names=['chrom', 'start', 'end', 'reads'], \
                            dtype={'chrom':'str', 'start':'int', 'end':'int', 'reads':'float'})

        if WHAT_TO_PLOT == 'tp':
            tp_file = "%s/%s_RNN-hybrid_true_positives.bed" % (DATA_DIR, sample)

            tp_df = pd.read_csv(tp_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                                names=['chrom', 'start', 'end'],
                                dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
            tp_df = tp_df[tp_df['chrom'] == chromosome]
            count += len(tp_df)

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
            fp_file = "%s/%s_RNN-hybrid_false_positives.bed" % (DATA_DIR, sample)

            fp_df = pd.read_csv(fp_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                                names=['chrom', 'start', 'end'],
                                dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
            fp_df = fp_df[fp_df['chrom'] == chromosome]
            count += len(fp_df)

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
            tn_file = "%s/%s_RNN-hybrid_true_negatives.bed" % (DATA_DIR, sample)

            tn_df = pd.read_csv(tn_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                                names=['chrom', 'start', 'end'],
                                dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
            tn_df = tn_df[tn_df['chrom'] == chromosome]
            count += len(tn_df)

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
            fn_file = "%s/%s_RNN-hybrid_false_negatives.bed" % (DATA_DIR, sample)

            fn_df = pd.read_csv(fn_file, header=None, sep="\t", na_filter=False, usecols=[0, 1, 2],
                                names=['chrom', 'start', 'end'],
                                dtype={ 'chrom':'str', 'start':'int', 'end  ':'int'})
            fn_df = fn_df[fn_df['chrom'] == chromosome]
            count += len(fn_df)

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
    return tp_pos_signal, tp_neg_signal, fp_pos_signal, fp_neg_signal, tn_pos_signal, tn_neg_signal, fn_pos_signal, fn_neg_signal, count
          



if __name__=='__main__':

    CHROMOSOMES = [ 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
                    'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    pool = multiprocessing.Pool(NR_THREADS)
    final_result = pool.map(get_metapeaks, CHROMOSOMES)
    pool.close()
    pool.join()
    final_result = np.array(final_result)
    aggregated = final_result.sum(axis=0)

    my_dpi = 300

    if WHAT_TO_PLOT == 'tp':
        np.save('%s/multi_nascent_TP_pos.npy' % DATA_DIR, aggregated[0])
        np.save('%s/multi_nascent_TP_neg.npy' % DATA_DIR, aggregated[1])
        plt.clf()
        fig = plt.figure(figsize=(1200/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[0], aggregated[8]), color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[1], aggregated[8]), color='red')
        plt.title('Nascent transcription signal meta-peak \n from true positives (n=%d)' % aggregated[8])
        plt.ylim((-0.12,0.12))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total normalized number of reads')
        plt.text(400, 0.08, "TP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
        plt.tight_layout()
        plt.savefig("%s/multi_TP_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'fp':
        np.save('%s/multi_nascent_FP_pos.npy' % DATA_DIR, aggregated[2])
        np.save('%s/multi_nascent_FP_neg.npy' % DATA_DIR, aggregated[3])
        plt.clf()
        fig = plt.figure(figsize=(1200/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[2], aggregated[8]), color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[3], aggregated[8]), color='red')
        plt.title('Nascent transcription signal meta-peak \n from false positives (n=%d)' % aggregated[8])
        plt.ylim((-0.12,0.12))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total normalized number of reads')
        plt.text(400, 0.08, "FP", horizontalalignment='left', size='xx-large', color='black', weight='bold')
        plt.tight_layout()
        plt.savefig("%s/multi_FP_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'tn':
        np.save('%s/multi_nascent_TN_pos.npy' % DATA_DIR, aggregated[4])
        np.save('%s/multi_nascent_TN_neg.npy' % DATA_DIR, aggregated[5])
        plt.clf()
        fig = plt.figure(figsize=(1200/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[4], aggregated[8]), color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[5], aggregated[8]), color='red')
        plt.title('Nascent transcription signal meta-peak \n from true negatives (n=%d)' % aggregated[8])
        plt.ylim((-0.12,0.12))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total normalized number of reads')
        plt.text(400, 0.08, "TN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
        plt.tight_layout()
        plt.savefig("%s/multi_TN_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'fn':
        np.save('%s/multi_nascent_FN_pos.npy' % DATA_DIR, aggregated[6])
        np.save('%s/multi_nascent_FN_neg.npy' % DATA_DIR, aggregated[7])
        plt.clf()
        fig = plt.figure(figsize=(1200/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[6], aggregated[8]), color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(aggregated[7], aggregated[8]), color='red')
        plt.title('Nascent transcription signal meta-peak \n from false negatives (n=%d)' % aggregated[8])
        plt.ylim((-0.12,0.12))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total normalized number of reads')
        plt.text(400, 0.08, "FN", horizontalalignment='left', size='xx-large', color='black', weight='bold')
        plt.tight_layout()
        plt.savefig("%s/multi_FN_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'all_p':
        tp_pos = np.load('%s/multi_nascent_TP_pos.npy' % DATA_DIR)
        tp_neg = np.load('%s/multi_nascent_TP_neg.npy' % DATA_DIR)
        fn_pos = np.load('%s/multi_nascent_FN_pos.npy' % DATA_DIR)
        fn_neg = np.load('%s/multi_nascent_FN_neg.npy' % DATA_DIR)

        tp_count = 0
        for sample in ['A549', 'GM12878', 'H1', 'HeLa', 'K562', 'LNCaP', 'MCF7', 'THP1']:
            tp_file = "%s/%s_RNN-hybrid_true_positives.bed" % (DATA_DIR, sample)
            with open(tp_file, 'r') as tp_fd:
                for _ in tp_fd.readlines():
                    tp_count += 1

        fn_count = 0
        for sample in ['A549', 'GM12878', 'H1', 'HeLa', 'K562', 'LNCaP', 'MCF7', 'THP1']:
            fn_file = "%s/%s_RNN-hybrid_false_negatives.bed" % (DATA_DIR, sample)
            with open(fn_file, 'r') as fn_fd:
                for _ in fn_fd.readlines():
                    fn_count += 1

        total_count = tp_count + fn_count
        allp_pos = tp_pos + fn_pos
        allp_neg = tp_neg + fn_neg

        mpl.rc('axes', edgecolor='green')
        mpl.rc('xtick', color='green')
        mpl.rc('ytick', color='green')
        plt.clf()
        fig = plt.figure(figsize=(1200/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(allp_pos, total_count), color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(allp_neg, total_count), color='red')
        plt.title('Nascent transcription signal meta-peak \n from all positives (n=%d)' % total_count)
        plt.ylim((-0.1,0.1))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total normalized number of reads')
        plt.tight_layout()
        plt.savefig("%s/multi_allP_nascent_metapeak.png" % DATA_DIR)

    elif WHAT_TO_PLOT == 'all_n':
        tn_pos = np.load('%s/multi_nascent_TN_pos.npy' % DATA_DIR)
        tn_neg = np.load('%s/multi_nascent_TN_neg.npy' % DATA_DIR)
        fp_pos = np.load('%s/multi_nascent_FP_pos.npy' % DATA_DIR)
        fp_neg = np.load('%s/multi_nascent_FP_neg.npy' % DATA_DIR)

        tn_count = 0
        for sample in ['A549', 'GM12878', 'H1', 'HeLa', 'K562', 'LNCaP', 'MCF7', 'THP1']:
            tn_file = "%s/%s_RNN-hybrid_true_negatives.bed" % (DATA_DIR, sample)
            with open(tn_file, 'r') as tn_fd:
                for _ in tn_fd.readlines():
                    tn_count += 1

        fp_count = 0
        for sample in ['A549', 'GM12878', 'H1', 'HeLa', 'K562', 'LNCaP', 'MCF7', 'THP1']:
            fp_file = "%s/%s_RNN-hybrid_false_positives.bed" % (DATA_DIR, sample)
            with open(fp_file, 'r') as fp_fd:
                for _ in fp_fd.readlines():
                    fp_count += 1

        total_count = tn_count + fp_count
        alln_pos = tn_pos + fp_pos
        alln_neg = tn_neg + fp_neg

        mpl.rc('axes', edgecolor='green')
        mpl.rc('xtick', color='green')
        mpl.rc('ytick', color='green')
        plt.clf()
        fig = plt.figure(figsize=(1200/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.zeros(2*EVALUATION_RADIUS), color='black')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(alln_pos, total_count), color='blue')
        plt.step(list(range(-1*EVALUATION_RADIUS,EVALUATION_RADIUS)), np.divide(alln_neg, total_count), color='red')
        plt.title('Nascent transcription signal meta-peak \n from all negatives (n=%d)' % total_count)
        plt.ylim((-0.1,0.1))
        plt.xlabel('Nucleotide position from midpoint')
        plt.ylabel('Total normalized number of reads')
        plt.tight_layout()
        plt.savefig("%s/multi_allN_nascent_metapeak.png" % DATA_DIR)

