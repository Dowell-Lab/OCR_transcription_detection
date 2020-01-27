import pandas as pd

data = pd.read_pickle('./combined_dataset_union_fstitchtfit.pkl')
data['prom_ovlp'] = data['prom_ovlp'].clip(upper=1)

print("For all OCRs:")
print("OCRs overlapping nascent transcription:")
print(data.groupby(['sample', 'ovlp_txn']).size())
print("OCRs overlapping TSSs:")
print(data.groupby(['sample', 'prom_ovlp']).size())

print("For the test set (no HCT116, chr12-chrY only):")
chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY']
split_index = 11 # ie. chr12
test_chroms = chromosomes[split_index:]
data = data[(data['chrom'].isin(test_chroms)) & (data['sample'] != 'HCT116')]

print("OCRs overlapping nascent transcription:")
print(data.groupby(['sample', 'ovlp_txn']).size())
print("OCRs overlapping TSSs:")
data['prom_ovlp'].clip(upper=1)
print(data.groupby(['sample', 'prom_ovlp']).size())



