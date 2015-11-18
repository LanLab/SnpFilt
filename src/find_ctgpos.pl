#!/usr/bin/env perl

die 
"Usage: ./find_ctgpos.pl <cigarxfile> <posfile> <out1-file> <out2-file>
Input file posfile should in csv format with each line in format: 
    <chrom>,<position>
Two output files are generated
    <out1-file> will have format <ctg>,<cpos> - this is used with zgrep so matching needs to be exact
    <out2-file> will have format <chrom><gpos>,<ctg>,<cpos>,<cdir> - will be used for merging\n"
if $#ARGV != 3;

$alnfile = $ARGV[0];
$posfile = $ARGV[1];
$out1file = $ARGV[2];
$out2file = $ARGV[3];

$MIN_QUAL = 20;
$MIN_MQ = 30;
$MIN_NREADS = 10;
$MIN_NREAD_DIR = 1;
$MIN_HETERO = 0.75;

$COL_CTG = 4;
$COL_CHROM = 0;
$COL_GSTART = 1;
$COL_GEND = 2;
$COL_LSTART = 5;
$COL_LEND = 6;
$COL_LSTRAND = 7;
$COL_CIGARX = 13;   # NOT 14!! 


# Read in alignment file (generated using lastz)
my @chrom = ();
my @ctgname = ();
my @lstart = ();
my @lend = ();
my @lstrand = ();
my @gstart = ();
my @gend = ();
my @cigarx = ();

open ALN, "<", $alnfile or die "Cannot open $alnfile for reading\n";
$line = <ALN>;  # skip header line
while ($line = <ALN>)
{
    @ln = split /\t/, $line;
    push @chrom, $ln[$COL_CHROM];
    push @ctgname, $ln[$COL_CTG];
    push @lstart, $ln[$COL_LSTART];
    push @lend, $ln[$COL_LEND];
    push @lstrand, $ln[$COL_LSTRAND];
    push @gstart, $ln[$COL_GSTART];
    push @gend, $ln[$COL_GEND];
    push @cigarx, $ln[$COL_CIGARX];
}
close ALN;
@ind = sort {$chrom[$a] cmp $chrom[$b] || $gstart[$a] <=> $gstart[$b]} (0..$#ctgname);

# Read in csv-file containing genomic positions of interest
# Sort and re-order arrays (but they are assumed to be unique)
@pchrom = ();
@ppos = ();
open POS, "<", $posfile or die "Cannot open $posfile\n";
while ($line = <POS>)
{
    next if substr($line,0,1) eq '#';

    chomp $line;
    @ln = split /,/, $line;
    push @pchrom, $ln[0];
    push @ppos, $ln[1];
}
close POS;
@g_ordered = sort {$pchrom[$a] cmp $pchrom[$b] || $ppos[$a] <=> $ppos[$b]} (0..$#ppos);
@pchrom = @pchrom[@g_ordered];
@ppos = @ppos[@g_ordered];
    
#foreach $ppi (0..$#pchrom)
#{
#    print "pchrom=$pchrom[$ppi], gpos=$ppos[$ppi]\n";
#}
#exit;

# Iterate through alignments for SNPs
$ppi = 0;

open OUT1, ">", $out1file or die "Cannot open $out1file for writing\n";
open OUT2, ">", $out2file or die "Cannot open $out2file for writing\n";
print OUT1 "ctg,cpos\n";
print OUT2 "chrom,gpos,ctg,cpos,cdir\n";
while ($ppi <= $#ppos)
{
foreach $i (0..$#ind)
{

    $tmp_cigarx = $cigarx[$ind[$i]];
    $boffset = 0;
    $goffset = 0;
    $loffset = 0;

    while ($ppi <= $#ppos && $ppos[$ppi] < $gstart[$ind[$i]])
    {
        $ppi++;
    }
    next if $ppos[$ppi] > $gend[$ind[$i]];

    while ($tmp_cigarx=~/(\d*)([MID])/g)
    {
        $slen = $1 + 0;
        $stype = uc($2);
        $slen = 1 if $slen == 0;
        $snp_gpos = $gstart[$ind[$i]] + $boffset + $goffset + $slen;
        if ($stype eq "M") 
        {
            while ($ppi <= $#ppos && $ppos[$ppi] <= $snp_gpos)
            {
                $poffset = $ppos[$ppi] - $gstart[$ind[$i]];
                $lpos = $lstart[$ind[$i]] + $loffset - $goffset + $poffset;
                $lpos = ($lend[$ind[$i]] + 1) - ($loffset - $goffset + $poffset) if $lstrand[$ind[$i]] eq '-';

                #print "$ppi $pchrom[$ppi] $ppos[$ppi] $ctgname[$ind[$i]] $lpos $lstrand[$ind[$i]]\n";
                print OUT1 "$ctgname[$ind[$i]],$lpos\n";
                print OUT2 "$pchrom[$ppi],$ppos[$ppi],$ctgname[$ind[$i]],$lpos,$lstrand[$ind[$i]]\n";
                $ppi++;
            }
        }
        $boffset += $slen if $stype eq 'M';
        $loffset += $slen if $stype eq 'I';
        $goffset += $slen if $stype eq 'D';
    }
}
}   


