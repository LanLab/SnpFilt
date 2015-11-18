#!/usr/bin/env perl

die "Usage: best_match.pl <cigarxfile>\n" if $#ARGV != 0;

$infile = $ARGV[0];

# Variables for the alignment
$MIN_ALN_LEN = 1200;
$MAX_OVERLAP = 100;
$MIN_IDPCT = 0.00;

# Format-specific constant
$COL_CTG = 4;
$COL_LSTART = 5;
$COL_LEND = 6;
$COL_SCORE = 8;
$COL_IDPCT = 10;

open CGX, "<", $infile or die "Cannot open $infile for reading\n";
@cgxdat = <CGX>;
close CGX;

%ctgct = ();
@score = (-1)x($#cgxdat+1);
@lstart = (-1)x($#cgxdat+1);
@lend = (-1)x($#cgxdat+1);

foreach $ii (1..$#cgxdat)
{
    @ln = split /\t/, $cgxdat[$ii];
    $lstart[$ii] = $ln[$COL_LSTART];
    $lend[$ii] = $ln[$COL_LEND];
    $idpct = $ln[$COL_IDPCT] + 0.0;

    next if ($lend[$ii] - $lstart[$ii] < $MIN_ALN_LEN);
    next if ($idpct < $MIN_IDPCT);
    $score[$ii] = $ln[$COL_SCORE];
    $ctgct{$ln[$COL_CTG]}++;
}

print $cgxdat[0];
foreach $ctg (sort keys %ctgct)
{
    # Skip if there are no valid alignments
    # These contigs will still be retained in merge_contigs.pl
    next if $ctgct{$ctg} == 0;

    @ind0 = grep {$score[$_] > 0 && $cgxdat[$_]=~/$ctg/} (1..$#cgxdat);
    @ind = sort {$lstart[$a] <=> $lstart[$b]} @ind0;
    @keep = (1)x($#ind+1);
    
    for ($i = 0; $i <= $#ind; $i++) 
    {
        next if $keep[$i] == 0;
        for ($j = $i+1; $j <= $#ind && ($lend[$ind[$i]] - $lstart[$ind[$j]]) > $MAX_OVERLAP; $j++)
        {
            if ($score[$ind[$j]] > $score[$ind[$i]])
            {
                $keep[$i] = 0;
                last;
            }
            else
            {
                $keep[$j] = 0;
            }
        }
    }
   
    foreach $i (0..$#ind)
    {
        print  $cgxdat[$ind[$i]] if $keep[$i];
    }
}
close OUT;






