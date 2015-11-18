#!/usr/bin/env perl

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
$COL_CIGARX = 14; 

die "Usage: get_seqvar_pos.pl <cigarxfile> <fasta1> <fasta2>\n" if $#ARGV != 2;
$alnfile = $ARGV[0];

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

# read in sequences
%seq1 = ();
%seq2 = ();

foreach $i (1,2)
{
    $seq_href = \%seq1 if $i == 1;
    $seq_href = \%seq2 if $i == 2;

    $seqfile = $ARGV[$i];
    open SEQ, "<", $seqfile or die "Cannot open $seqfile\n";
    $tmpseq = "";
    $tmpname = "";
    while ($line = <SEQ>)
    {
        chomp $line;
        if (substr($line, 0, 1) eq '>' || eof(SEQ))
        {
            $seq_href->{$tmpname} = $tmpseq if (length($tmpseq) > 0);
            $tmpseq = "";
            $tmpname = $line;
            $tmpname =~s/^>\s*//;
            $tmpname =~s/\s+.*$//;
        }
        else
        {
            $tmpseq .= $line;
        }
    }
    $seq_href->{$tmpname} = $tmpseq;
    close SEQ;
}

# Iterate through alignments for SNPs
print "chrom,gpos,ref,alt,stype,slen,ctg,cpos,cstrand\n";
foreach $i (0..$#ind)
{
    $tmp_cigarx = $cigarx[$ind[$i]];
    $boffset = 0;
    $goffset = 0;
    $loffset = 0;
    while ($tmp_cigarx=~/(\d*)([=XID])/g)
    {
        $slen = $1 + 0;
        $stype = uc($2);
        $slen = 1 if $slen == 0;
        if ($stype ne '=')
        {
            $snp_gpos = $gstart[$ind[$i]] + $boffset + $goffset + 1;
            $snp_lpos = $lstart[$ind[$i]] + $boffset + $loffset + 1;
            $snp_lpos = $lend[$ind[$i]] - ($boffset + $loffset + $slen) + 1 if $lstrand[$ind[$i]] eq '-';
                
            #$ref = substr($seq1{$chrom[$ind[$i]]}, $snp_gpos-1, $slen) if ($stype ne 'I');
            #$alt = substr($seq2{$ctgname[$ind[$i]]}, $snp_lpos-1, $slen) if ($stype ne 'D');

            $snp_ext = 1; 
            $snp_ext = $slen if $stype eq 'X';
            $snp_len2 = $slen;
            $snp_len2 = 1 if $stype eq 'X'; 
            $cdir = 1;
            $cdir = -1 if $lstrand[$ind[$i]] eq '-';

            for ($snp_off = 0; $snp_off < $snp_ext; $snp_off++)
            {
                $snp_gpos2 = $snp_gpos + $snp_off;
                $snp_lpos2 = $snp_lpos + ($cdir * $snp_off);
                $ref = "-";
                $alt = "-";
                $ref = substr($seq1{$chrom[$ind[$i]]}, $snp_gpos2-1, $snp_len2) if $stype ne 'I';
                $alt = substr($seq2{$ctgname[$ind[$i]]}, $snp_lpos2-1, $snp_len2) if $stype ne 'D';
                
                if ($lstrand[$ind[$i]] eq '-')
                {
                    $alt = reverse($alt);
                    $alt =~tr/ATCG/TAGC/;
                }
                print "$chrom[$ind[$i]],$snp_gpos2,$ref,$alt,$stype,$slen,$ctgname[$ind[$i]],$snp_lpos2,$lstrand[$ind[$i]]\n";
            } 
        }
        
        $boffset += $slen if $stype eq '=';
        $boffset += $slen if $stype eq 'X';
        $loffset += $slen if $stype eq 'I';
        $goffset += $slen if $stype eq 'D';
    }
}
   

