#!/bin/sh

cd ${OUTDIR}

MAPPED=`cat ./${1}/qc/mapstats/${2}.millionsmapped`
MILMAPPED=`echo "scale=8 ; ${MAPPED} / 1000000" | bc`
awk -v millions="${MILMAPPED}" 'BEGIN {FS = "\t";OFS = "\t"} {$4 = $4 / millions}1' ./${1}/mapped/bedgraphs/${2}.bedGraph > ${OUTDIR}/${2}.rcc.BedGraph
