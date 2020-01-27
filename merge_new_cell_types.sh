#!/bin/sh


cd ${BEDGRAPHS}

# new cell types:
declare -a arr=( "SRR1810067" "SRR1810068" "SRR1810069" "SRR1810070" "SRR1810071" "SRR1810072" "SRR1648886" "SRR1648890" "SRR1648896" "SRR1648897" "SRR1648909" )

# split the normalized bedGraphs into positive and negative strand files
for S in "${arr[@]}"
do
    grep -v "-" ${S}.rcc.bedGraph > ${S}.pos.rcc.bedGraph
    grep "-" ${S}.rcc.bedGraph > ${S}.neg.rcc.bedGraph
done


# New cell types
bedtools unionbedg -i SRR1648886.pos.rcc.bedGraph \
                      SRR1648890.pos.rcc.bedGraph \
                      SRR1648896.pos.rcc.bedGraph \
                      SRR1648897.pos.rcc.bedGraph \
                      SRR1648909.pos.rcc.bedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > lncap.pos.bedGraph

bedtools unionbedg -i SRR1648886.neg.rcc.bedGraph \
                      SRR1648890.neg.rcc.bedGraph \
                      SRR1648896.neg.rcc.bedGraph \
                      SRR1648897.pos.rcc.bedGraph \
                      SRR1648909.pos.rcc.bedGraph \
  | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' \
  > lncap.neg.bedGraph

sed 's/-//g' lncap.neg.bedGraph > lncap.fakeneg.bedGraph

bedtools unionbedg -i lncap.pos.bedGraph lncap.fakeneg.bedGraph | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > lncap.nascent.bedGraph


