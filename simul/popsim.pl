#!/usr/bin/perl -w
use strict;

my $USAGE = "popsim.pl num_founders num_samples var_per_founder var_per_sample > pop.var\n";

die $USAGE if (scalar @ARGV != 4);

my $NUMFOUNDERS = shift @ARGV;
my $NUMSAMPLES  = shift @ARGV;
my $VARFOUNDER  = shift @ARGV;
my $VARSAMPLE   = shift @ARGV;

my $CHRLEN = 100000000;

print STDERR "Simulating $NUMFOUNDERS founders with $NUMSAMPLES samples each. Founders have $VARFOUNDER variants, samples have $VARSAMPLE extra\n";

my %popvar;

my @samples;

for (my $founder = 0; $founder < $NUMFOUNDERS; $founder++)
{
  print STDERR " simulating founder $founder\n";
  my @foundervars;

  for (my $i = 0; $i < $VARFOUNDER; $i++)
  {
    my $pos = int(rand($CHRLEN));
    push @foundervars, $pos;
  }

  for (my $sample = 0; $sample < $NUMSAMPLES; $sample++)
  {
    # Load the founder variants for this sample
    my $sampleid = "f${founder}s${sample}";
    push @samples, $sampleid;

    foreach my $pos (@foundervars)
    {
      $popvar{$pos}->{$sampleid} = 1;
    }

    # now pick the extra variants
    for (my $i = 0; $i < $VARSAMPLE; $i++)
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
print "##reference=${NUMFOUNDERS}founders:${NUMSAMPLES}samples:${VARFOUNDER}foundervar:${VARSAMPLE}samplevar\n";
print "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach my $sampleid (@samples) { print "\t$sampleid"; }
print "\n";

foreach my $pos (sort {$a <=> $b} keys %popvar)
{
  print "chr1\t$pos\t.\ta\t<INS>\t60\tPASS\tNS=$totalsamples\tGT";

  foreach my $sampleid (@samples)
  {
    my $gt = "0/0";
    if (exists ($popvar{$pos}->{$sampleid})) { $gt = "1/1"; }

    print "\t$gt";
  }
  print "\n";
}
