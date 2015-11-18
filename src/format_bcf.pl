#!/usr/bin/env perl

die "Usage: format_bcf.pl <bcffile> <out.gz>\n" if $#ARGV != 1;
$bcffile = $ARGV[0];
$outfile = $ARGV[1];

die "error: environment varable BCFTOOLS not set\n" if (!exists($ENV{'BCFTOOLS'}));

open OUT, '|-', "gzip >$outfile" or die "Cannot open gzip $outfile\n";
print OUT "ctg,cpos,ref,qual,mq,dpFR,dpRR,dptot\n";
open BCF, "$ENV{'BCFTOOLS'} call -c -V indels $bcffile |" or die "Cannot open $bcffile\n";
while ($line = <BCF>)
{
    next if substr($line, 0, 1) eq '#';

    @ln = split /\t/, $line;
    $ln[7] =~/DP4=(\d+),(\d+),(\d+),(\d+)/;
    @dp = ($1,$2,$3,$4);
    $dptot=$dp[0]+$dp[1]+$dp[2]+$dp[3];
    $ln[7] =~/MQ=(\d+)/; $mq = $1;
    print OUT "$ln[0],$ln[1],$ln[3],$ln[5],$mq,$dp[0],$dp[1],$dptot\n";
}
close OUT;
close BCF;

