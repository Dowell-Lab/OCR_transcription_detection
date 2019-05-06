#!/bin/bash

cd /scratch/Shares/dowell/processed_atacseq/macs2


cp SRR1822165_peaks_clean.broadPeak gm12878_peaks_union.broadPeak
cat SRR1822166_peaks_clean.broadPeak >> gm12878_peaks_union.broadPeak
cat SRR1822167_peaks_clean.broadPeak >> gm12878_peaks_union.broadPeak
cat SRR1822168_peaks_clean.broadPeak >> gm12878_peaks_union.broadPeak
bedtools sort -i gm12878_peaks_union.broadPeak > gm12878_peaks_union.sorted.broadPeak
bedtools merge -i gm12878_peaks_union.sorted.broadPeak > gm12878_peaks_union.broadPeak
rm gm12878_peaks_union.sorted.broadPeak

cp SRR5007258_peaks_clean.broadPeak h1_peaks_union.broadPeak
cat SRR5007259_peaks_clean.broadPeak >> h1_peaks_union.broadPeak
bedtools sort -i h1_peaks_union.broadPeak > h1_peaks_union.sorted.broadPeak
bedtools merge -i h1_peaks_union.sorted.broadPeak > h1_peaks_union.broadPeak
rm h1_peaks_union.sorted.broadPeak

cp SRR5876158_peaks_clean.broadPeak hct116_peaks_union.broadPeak
cat SRR5876159_peaks_clean.broadPeak >> hct116_peaks_union.broadPeak
bedtools sort -i hct116_peaks_union.broadPeak > hct116_peaks_union.sorted.broadPeak
bedtools merge -i hct116_peaks_union.sorted.broadPeak > hct116_peaks_union.broadPeak
rm hct116_peaks_union.sorted.broadPeak

cp SRR5128074_peaks_clean.broadPeak k562_peaks_union.broadPeak

echo "Processing CD4+"
cp SRR891275_peaks_clean.broadPeak cd4pos_peaks_union.broadPeak
cat SRR891276_peaks_clean.broadPeak >> cd4pos_peaks_union.broadPeak
bedtools sort -i cd4pos_peaks_union.broadPeak > cd4pos_peaks_union.sorted.broadPeak
bedtools merge -i cd4pos_peaks_union.sorted.broadPeak > cd4pos_peaks_union.broadPeak
rm cd4pos_peaks_union.sorted.broadPeak

echo "Processing Jurkat"
cp SRR5063984_peaks_clean.broadPeak jurkat_peaks_union.broadPeak
cat SRR5063985_peaks_clean.broadPeak >> jurkat_peaks_union.broadPeak
cat SRR5063986_peaks_clean.broadPeak >> jurkat_peaks_union.broadPeak
bedtools sort -i jurkat_peaks_union.broadPeak > jurkat_peaks_union.sorted.broadPeak
bedtools merge -i jurkat_peaks_union.sorted.broadPeak > jurkat_peaks_union.broadPeak
rm jurkat_peaks_union.sorted.broadPeak

cp SRR3622817_peaks_clean.broadPeak lncap_peaks_union.broadPeak
cat SRR3622818_peaks_clean.broadPeak >> lncap_peaks_union.broadPeak
cat SRR3622819_peaks_clean.broadPeak >> lncap_peaks_union.broadPeak
bedtools sort -i lncap_peaks_union.broadPeak > lncap_peaks_union.sorted.broadPeak
bedtools merge -i lncap_peaks_union.sorted.broadPeak > lncap_peaks_union.broadPeak
rm lncap_peaks_union.sorted.broadPeak


# Merge bedGraph files, set the depth column to the average of all samples for that region
#cd /scratch/Shares/dowell/processed_atacseq/bedtools/

bedtools unionbedg -i SRR1822165.sorted2.mp.BedGraph \
                      SRR1822167.sorted2.mp.BedGraph \
                      SRR1822167.sorted2.mp.BedGraph \
                      SRR1822168.sorted2.mp.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > gm12878.sorted2.mp.BedGraph

bedtools unionbedg -i SRR5007258.sorted2.mp.BedGraph \
                      SRR5007259.sorted2.mp.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > h1.sorted2.mp.BedGraph

bedtools unionbedg -i SRR5876158.sorted2.mp.BedGraph \
                      SRR5876159.sorted2.mp.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > hct116.sorted2.mp.BedGraph

cp SRR5128074.sorted2.mp.BedGraph k562.sorted2.mp.BedGraph

bedtools unionbedg -i SRR891275.sorted2.mp.BedGraph \
                      SRR891276.sorted2.mp.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > cd4pos.sorted2.mp.BedGraph

bedtools unionbedg -i SRR5063984.sorted2.mp.BedGraph \
                      SRR5063985.sorted2.mp.BedGraph \
                      SRR5063986.sorted2.mp.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > jurkat.sorted2.mp.BedGraph

bedtools unionbedg -i SRR3622817.sorted.mp.BedGraph \
                      SRR3622818.sorted2.mp.BedGraph \
                      SRR3622819.sorted2.mp.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > lncap.sorted2.mp.BedGraph

