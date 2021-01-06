#!/bin/sh

if [ $# -ne 3 ]
then
  echo "USAGE: SVCollector.sh samples.vcf numtoplot workdir"
  exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SVC=$DIR/Debug/SVCollector 
PLOT=$DIR/vis/SVCollectorPlot.R

VCF=$1
NUMTOPLOT=$2
OUTDIR=$3

mkdir -p $OUTDIR
VCFB=`basename $VCF`

RANDOMTRIALS=10

echo "Analyzing $VCF"

echo " computing greedy selection (disable min allele freq)"
$SVC greedy $VCF -1 $NUMTOPLOT 0 $OUTDIR/$VCFB.greedy >& $OUTDIR/$VCFB.greedy.log

echo " computing topN selection (disable min allele freq)"
$SVC topN $VCF $NUMTOPLOT 0 $OUTDIR/$VCFB.topN >& $OUTDIR/$VCFB.topN.log

echo " computing random selection over $RANDOMTRIALS trials (disable min allele freq) "
(seq $RANDOMTRIALS | parallel -t $SVC random $VCF $NUMTOPLOT 0 $OUTDIR/$VCFB.random.{}) >& $OUTDIR/$VCFB.randomlog

echo " cleaning _tmp files"
rm -f $OUTDIR/*_tmp

echo " plotting results"
$PLOT $OUTDIR/$VCFB.greedy $OUTDIR/$VCFB.topN $OUTDIR/$VCFB.random $OUTDIR/$VCFB.png 
