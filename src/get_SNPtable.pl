#/usr/bin/env perl
use Getopt::Long;

%REV_BASE = ('C'=>'G', 'G'=>'C', 'A'=>'T', 'T'=>'A');

die
"Usage: get_SNPtable.pl <fstem> <topdir>
This scripts expects the files aln_<fstem>.cigarx, snpfilt_results_<fstem>.csv and flt_contigs.gz
in each subdirectory inside <topdir>\n" if $#ARGV != 1;

die "Error: environment variable SCRIPTDIR undefined\n" if (!exists($ENV{'SCRIPTDIR'}));

# File names
$fstem  = $ARGV[0];
$topdir = $ARGV[1];
@dirs = `ls -d $topdir/*/`;
die "No directories inside $topdir" if $#dirs==0;
$allpos = "all_snppos.csv";

system ("grep \'^0,\' $topdir/*/snpfilt_results\_$fstem.csv | grep \',X\' | cut -d, -f 2,3 | sort | uniq > $allpos");
foreach $i (0..$#dirs)
{
    $dd = $dirs[$i];
    chomp $dd;
    $alnfile = "$dd/aln_$fstem.cigarx";
    $fltfile = "$dd/flt_contigs.gz";
    $resfile = "$dd/snpfilt_results\_$fstem.csv";
    die "snpfilt.pl: get_SNPtable.pl: files $alnfile, $fltfile or $resfile emissing in $dd\n"
          unless (-e $alnfile && -e $fltfile && -e $resfile);

    $dirs[$i] =~s/^[^\w]*//;
    $dirs[$i] =~s/[^\w]*$//;
    $dd = $dirs[$i];
    system ("$ENV{'SCRIPTDIR'}/find_ctgpos.pl $alnfile $allpos cpos_$dd.tmp cgpos_$dd.tmp\n");
    system("zgrep -w -f cpos_$dd.tmp $fltfile > cposF_$dd.tmp");
}

# Merge all files into one table
open POS, "<", $allpos or die "Cannot open $allpos\n";
@chrom = ();
@gpos = ();
while ($line = <POS>)
{
    chomp $line;
    @ln = split /,/, $line;
    push @chrom, $ln[0];
    push @gpos, $ln[1];
}
close POS;
@ord = sort {@chrom[$a] cmp $chrom[$b] || $gpos[$a] <=> $gpos[$b]} (0..$#gpos);
@chrom = @chrom[@ord];
@gpos = @gpos[@ord];

print ",", join(",", @chrom), "\n";
print ",", join(",", @gpos), "\n";

foreach $dd (@dirs)
{
    print "$dd";
    open INFO, "<", "cgpos_$dd.tmp" or die "Cannot open cgpos_$dd.tmp\n";
    @info = <INFO>;
    chomp @info;
    close INFO;
    foreach $i (0..$#gpos)
    {
        $seq = 'N';
        $gpos_str = quote("$chrom[$i],$gpos[$i],");
        @cpos_match = grep {$_=~/$gpos_str/} @info;
        @ln = split /,/, $cpos_match[0];
        $cpos_str = "$ln[2],$ln[3]";
        $cdir = $ln[4];
        if (length($cpos_str) > 2)
        {
            $cmd2 = "grep \'$cpos_str,\' cposF_$dd.tmp";
            $cft_line = `$cmd2`;
            chomp $cft_line;
        
            @cflt = split /,/, $cft_line;
            if ($cflt[2] == 0)
            {
                $seq = $cflt[3];
                $seq = $REV_BASE{$seq} if $cdir eq '-';
            }
        }
        print ",$seq"; 
    }
    print "\n";
}
    


 

