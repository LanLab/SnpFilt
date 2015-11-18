#!/usr/bin/env perl
#
die "Usage: convert_fasta.pl <flt-gzipped file>\n" if $#ARGV != 0;
$infile = $ARGV[0];

open IN, "gunzip -c $infile |" or die "Cannot open gzipped file $infile";
$line = <IN>;   # header line

$oldctg = "";
while ($line = <IN>)
{
    chomp $line;
    @ln = split /,/, $line;

    if ($ln[0] ne $oldctg)
    {
        print "\n" unless length($oldctg)==0;
        print ">$ln[0]\n";
        for ($i = 1; $i < $ln[1]; $i++)
        {
            print "N";
        }
    }
    else
    {
        for ($i = $oldpos + 1; $i < $ln[1]; $i++)
        {
            print "N";
        }
        $base = $ln[3];
        $base = "N" if $ln[2] != 0;
        print "$base";
    }
    $oldctg = $ln[0];
    $oldpos = $ln[1];
}
close IN;
    
