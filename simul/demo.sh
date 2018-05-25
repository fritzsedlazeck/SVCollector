#!/bin/sh

echo "Simulating and analyzing simple population"

## Create a simple population with 10 founder genomes
(./popsim.pl 10 10 0 1000 0 100 0 > simple10.vcf) >& simple10.vcf.log

## Now analyze the population focusing on the first 50 samples
../SVCollector.sh simple10.vcf 50 simple10

echo
echo

echo "Simulating and analyzing complex population"

## Create a complicated population with 10 founder genomes with an average of 10 samples per genome
(./popsim.pl 10 10 .5 500 .5 500 .5 > complex10.vcf) >& complex10.vcf.log

## Now analyze the population focusing on the first 50 samples
../SVCollector.sh complex10.vcf 50 complex10
