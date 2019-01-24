#!/usr/bin/env perl
use strict;
use Math::Random;

my $USAGE = "popsim.pl num_founders sample_per_founder std_samples var_per_founder std_per_founder var_per_sample std_per_sample > pop.var\n";

die $USAGE if (scalar @ARGV != 7);

my $NUMFOUNDERS = shift @ARGV;
my $NUMSAMPLES  = shift @ARGV;
my $STDSAMPLES  = shift @ARGV;

my $VARFOUNDER  = shift @ARGV;
my $STDFOUNDER  = shift @ARGV;

my $VARSAMPLE   = shift @ARGV;
my $STDSAMPLE   = shift @ARGV;


my $CHRLEN = 100000000;

print STDERR "Simulating $NUMFOUNDERS founders with $NUMSAMPLES +/- $STDSAMPLES% samples each\n";
print STDERR "Founders have $VARFOUNDER +/- $STDFOUNDER% variants, samples have $VARSAMPLE +/- $STDSAMPLE% extra\n";


my %popvar;

my @samples;

for (my $founder = 0; $founder < $NUMFOUNDERS; $founder++)
{
  my @foundervars;

  my $thisfounder = int(random_normal(1, $VARFOUNDER, $VARFOUNDER*$STDFOUNDER));
  if ($thisfounder < 1) { $thisfounder = 1; }

  for (my $i = 0; $i < $thisfounder; $i++)
  {
    my $pos = int(rand($CHRLEN));
    push @foundervars, $pos;
  }

  my $thispopulation = int(random_normal(1, $NUMSAMPLES, $NUMSAMPLES * $STDSAMPLES));
  if ($thispopulation < 1) { $thispopulation = 1; }

  print STDERR " simulating founder $founder with $thispopulation samples with $thisfounder variants\n";

  for (my $sample = 0; $sample < $thispopulation; $sample++)
  {
    # Load the founder variants for this sample
    my $sampleid = "f${founder}s${sample}";
    push @samples, $sampleid;

    foreach my $pos (@foundervars)
    {
      $popvar{$pos}->{$sampleid} = 1;
    }

    # now pick the extra variants
    my $thissample = int(random_normal(1, $VARSAMPLE, $VARSAMPLE*$STDSAMPLE));
    if ($thissample < 0) { $thissample = 0; }
    print STDERR "   $sampleid: $thisfounder + $thissample extra variants\n";

    for (my $i = 0; $i < $thissample; $i++)
    {
      my $pos = int(rand($CHRLEN));
      $popvar{$pos}->{$sampleid} = 1;
    }
  }
}

my $totalsamples = scalar @samples;
my ($day, $month, $year) = (localtime)[3,4,5];


print "##fileformat=VCFv4.0\n";
print "##fileDate=$year$month$day\n";
print "##reference=${NUMFOUNDERS}founders:${NUMSAMPLES}samples:${STDSAMPLES}std_samples:${VARFOUNDER}foundervar:${STDFOUNDER}:${VARSAMPLE}samplevar:${STDSAMPLE}stdsample\n";
print "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach my $sampleid (@samples) { print "\t$sampleid"; }
print "\n";

my $totalsv = 0;

foreach my $pos (sort {$a <=> $b} keys %popvar)
{
  print "chr1\t$pos\t.\ta\t<INS>\t60\tPASS\tNS=$totalsamples\tGT";
  $totalsv++;

  foreach my $sampleid (@samples)
  {
    my $gt = "0/0";
    if (exists ($popvar{$pos}->{$sampleid})) { $gt = "1/1"; }

    print "\t$gt";
  }
  print "\n";
}

print STDERR "Simulated $totalsamples samples with $totalsv distinct variants\n";
