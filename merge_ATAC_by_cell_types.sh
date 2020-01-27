#!/bin/sh


cd ${BEDGRAPHSDIR}

bedtools unionbedg -i SRR7140571.rcc.BedGraph \
                      SRR7140572.rcc.BedGraph \
                      SRR7140573.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > A549.sorted2.mp.BedGraph

bedtools unionbedg -i SRR1822165.rcc.BedGraph \
                      SRR1822167.rcc.BedGraph \
                      SRR1822167.rcc.BedGraph \
                      SRR1822168.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > GM12878.sorted2.mp.BedGraph

bedtools unionbedg -i SRR5007258.rcc.BedGraph \
                      SRR5007259.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > H1.sorted2.mp.BedGraph

bedtools unionbedg -i SRR5876158.rcc.BedGraph \
                      SRR5876159.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > HCT116.sorted2.mp.BedGraph

bedtools unionbedg -i SRR6216226.rcc.BedGraph \
                      SRR6216227.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > HeLa.sorted2.mp.BedGraph

bedtools unionbedg -i SRR3622817.rcc.BedGraph \
                      SRR3622818.rcc.BedGraph \
                      SRR3622819.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > LNCaP.sorted2.mp.BedGraph

bedtools unionbedg -i SRX6443488.rcc.BedGraph \
                      SRX6443489.rcc.BedGraph \
                      SRX6443490.rcc.BedGraph \
                      SRX6443491.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > MCF7.sorted2.mp.BedGraph

bedtools unionbedg -i SRR8932925.rcc.BedGraph \
                      SRR8932927.rcc.BedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > THP1.sorted2.mp.BedGraph

