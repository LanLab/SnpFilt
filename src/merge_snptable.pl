#!/usr/bin/env perl
#
die "Usage: merge_snptable.pl <raw-snpfile> <flt-snpfile>\n" if $#ARGV != 1;

$rawfile = $ARGV[0];
$fltfile = $ARGV[1];

%snpflt = {};
open FLT, "<", $fltfile or die "Cannot open $fltfile\n";
<FLT>;
while ($line = <FLT>)
{
    chomp $line;
    @ln = split /,/, $line;
    $snpflt{$ln[0].",".$ln[1]} = $ln[2];
}
close FLT;

open RAW, "<", $rawfile or die "Cannot open $rawfile\n";
$line = <RAW>;

print
"# The 'flt' field describes the reliability of each site
# All values greater than 0 indicate a potential source of error
#
# 0: base is reliable
# 1: base has excessively high coverage
# 2: base is in neighbourhood with low mapping quality (MQ < 58)
# 3: base is in neighbourhood with low coverage (dpTot < 20 or dpFR < 1 or dpRR < 1)
# 4: base is in neighbourhood with low forward coverage (dpFR < 10)
# 5: base has high heterogeneity (dpNonRef > 0.3*dpTot)
# 6: base is in neighbourhood with low base quality
# 9: base does not occur in bcf file (coverage=0, check BCF files if frequent)
# \n";

print "flt,$line";
while ($line = <RAW>)
{
    @ln = split /,/, $line;
    $key = $ln[6].",".$ln[7];
    $flt = $snpflt{$key};
    $flt = 9 if length($flt) == 0;
    print "$flt,$line";
}
close RAW


