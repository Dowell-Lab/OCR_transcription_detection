import pandas as pd
import numpy as np
import multiprocessing as mp

NR_THREADS = 2
DATA_DIR = '.'
nascent_dir = '/home/ignaciot/dowell_lab/rnn_paper_backup/nascent_bg'
nascent_pos_suffix = '.pos.merged.clean.sorted.bedGraph'
nascent_neg_suffix = '.neg.merged.clean.sorted.bedGraph'

data = pd.read_pickle('%s/combined_dataset_union_fstitchtfit.pkl' % DATA_DIR)

last_cell_type = ''
pos_df = None
neg_df = None

EVALUATION_RADIUS = 500
chroms_used = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 
               'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 
               'chr8', 'chr9', 'chrX', 'chrY']

def nascent_coverage_for_cell_type(cell_type):
    nascent_coverage = []
    pos_file = "%s/%s%s" % (nascent_dir, cell_type, nascent_pos_suffix)
    neg_file = "%s/%s%s" % (nascent_dir, cell_type, nascent_neg_suffix)
    cell_df = data[(data['sample'] == cell_type)].sort_values(["chrom", "start"])
    nr_ocrs = len(cell_df)
    increment = int(nr_ocrs/10)
    print("Processing %s cells (%d OCRs)" % (cell_type, nr_ocrs))
    with open(pos_file, 'r') as pos_fd:
        with open(neg_file, 'r') as neg_fd:
            pos_line = pos_fd.readline()
            pos_chunks = pos_line.split('\t')
            neg_line = neg_fd.readline()
            neg_chunks = neg_line.split('\t')
            last_pos_chunks = ['none', '0', '0', '0\n']
            last_neg_chunks = ['none', '0', '0', '0\n']
            ocr_count = 0
            for ocr in cell_df.itertuples():
                midpoint = ocr.start + int((ocr.end - ocr.start)/2)
                window_start = midpoint - EVALUATION_RADIUS
                window_end = midpoint + EVALUATION_RADIUS
                ocr_size = 2*EVALUATION_RADIUS
                peak_window_reads = np.zeros(ocr_size)

                # advance the index until the beginning of that peak (or the next region past it)
                # TODO: figure out an efficient, smarter way to do this (it could only be an
                #       issue if there is NO coverage for an entire chromosome)
                while pos_line and chroms_used.index(ocr.chrom) > chroms_used.index(pos_chunks[0]):
                    pos_line = pos_fd.readline()
                    last_pos_chunks = pos_chunks
                    pos_chunks = pos_line.split('\t')
                while pos_line and ocr.chrom == pos_chunks[0] and window_start > int(pos_chunks[2]):
                    pos_line = pos_fd.readline()
                    last_pos_chunks = pos_chunks
                    pos_chunks = pos_line.split('\t')

                while neg_line and chroms_used.index(ocr.chrom) > chroms_used.index(neg_chunks[0]):
                    neg_line = neg_fd.readline()
                    last_neg_chunks = neg_chunks
                    neg_chunks = neg_line.split('\t')
                while neg_line and ocr.chrom == neg_chunks[0] and window_start > int(neg_chunks[2]):
                    neg_line = neg_fd.readline()
                    last_neg_chunks = neg_chunks
                    neg_chunks = neg_line.split('\t')

                while pos_line and ocr.chrom == pos_chunks[0] and int(pos_chunks[1]) < window_end:
                    feat_start = int(pos_chunks[1]) - window_start + 1
                    feat_end = int(pos_chunks[2]) - window_start
                    if feat_end >= ocr_size:
                        feat_end = ocr_size - 1
                    for position in range(int(feat_start), int(feat_end)):
                        peak_window_reads[position] = float(pos_chunks[3][:-1])
                    pos_line = pos_fd.readline()
                    last_pos_chunks = pos_chunks
                    pos_chunks = pos_line.split('\t')
 
                while neg_line and ocr.chrom == neg_chunks[0] and int(neg_chunks[1]) < window_end:
                    feat_start = int(neg_chunks[1]) - window_start + 1
                    feat_end = int(neg_chunks[2]) - window_start
                    if feat_end >= ocr_size:
                        feat_end = ocr_size - 1
                    for position in range(int(feat_start), int(feat_end)):
                        peak_window_reads[position] += abs(float(neg_chunks[3][:-1]))
                    neg_line = neg_fd.readline()
                    last_neg_chunks = neg_chunks
                    neg_chunks = neg_line.split('\t')

                # Before calling this a zero-nascent coverage window, check if the last positive or 
                # negative entry in the bedGraph was already scanned by a nearby upstream OCR
                if np.mean(peak_window_reads) == 0 and ocr.ovlp_txn == 1:
                    if last_pos_chunks[2] != 0 and ocr.chrom == last_pos_chunks[0] and int(last_pos_chunks[2]) > window_start:
                        feat_start = int(last_pos_chunks[1]) - window_start + 1
                        feat_end = int(last_pos_chunks[2]) - window_start
                        if feat_end >= ocr_size:
                            feat_end = ocr_size - 1
                        for position in range(int(feat_start), int(feat_end)):
                            peak_window_reads[position] = float(last_pos_chunks[3][:-1])

                    if last_pos_chunks[2] != 0 and ocr.chrom == last_neg_chunks[0] and int(last_neg_chunks[2]) > window_start:
                        feat_start = int(last_neg_chunks[1]) - window_start + 1
                        feat_end = int(last_neg_chunks[2]) - window_start
                        if feat_end >= ocr_size:
                            feat_end = ocr_size - 1
                        for position in range(int(feat_start), int(feat_end)):
                            peak_window_reads[position] += abs(float(last_neg_chunks[3][:-1]))

                    if np.mean(peak_window_reads) == 0:
                        print("Ran into an OCR that shouldn't be zero-coverage on nascent, trying the hard way...")
                        print(ocr.sample, ocr.chrom, ocr.start, ocr.end)
                        chrom_data_pos = pd.read_csv("%s/%s.%s.pos.merged.clean.sorted.bedGraph" % (nascent_dir, ocr.sample, ocr.chrom), sep="\t", na_filter=False, usecols=[1, 2, 3], names=['start', 'end', 'coverage'], dtype={'start':'int', 'end':'int', 'coverage':'float'})
                        chrom_data_pos = chrom_data_pos[(chrom_data_pos['start'] >= window_start) & (chrom_data_pos['end'] <= window_end)]
                        for coverage_chunk in chrom_data_pos.itertuples():
                            feat_start = coverage_chunk.start - window_start + 1
                            feat_end = coverage_chunk.end - window_start
                            if feat_end >= ocr_size:
                                feat_end = ocr_size - 1
                            for position in range(int(feat_start), int(feat_end)):
                                peak_window_reads[position] += abs(coverage_chunk.coverage)
                        del chrom_data_pos
                        chrom_data_neg = pd.read_csv("%s/%s.%s.neg.merged.clean.sorted.bedGraph" % (nascent_dir, ocr.sample, ocr.chrom), sep="\t", na_filter=False, usecols=[1, 2, 3], names=['start', 'end', 'coverage'], dtype={'start':'int', 'end':'int', 'coverage':'float'})
                        chrom_data_neg = chrom_data_neg[(chrom_data_neg['start'] >= window_start) & (chrom_data_neg['end'] <= window_end)]
                        for coverage_chunk in chrom_data_neg.itertuples():
                            feat_start = coverage_chunk.start - window_start + 1
                            feat_end = coverage_chunk.end - window_start
                            if feat_end >= ocr_size:
                                feat_end = ocr_size - 1
                            for position in range(int(feat_start), int(feat_end)):
                                peak_window_reads[position] += abs(coverage_chunk.coverage)
                        del chrom_data_neg
                        print(np.mean(peak_window_reads))

                nascent_coverage.append(np.mean(peak_window_reads))
                ocr_count += 1
                if ocr_count % increment == 0:
                    print("Progress of %s cells: %.1f%%" % (cell_type, (float(ocr_count)*100)/nr_ocrs))
    cell_df['mean_nr_nascent_reads'] = nascent_coverage
    return(cell_df)

pool = mp.Pool(NR_THREADS)
results = pool.map(nascent_coverage_for_cell_type, data['sample'].unique())
pool.close()
pool.join()

results_df = pd.concat(results)
results_df.to_pickle('%s/combined_dataset_union_fstitchtfit_with_nascent.pkl' % DATA_DIR)

