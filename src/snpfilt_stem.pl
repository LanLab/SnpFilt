
use Getopt::Long; 

$lastzscorefile="$SCRIPTDIR/lt2_dt104.scores";
$lastzargs="--nochain --gfextend --gapped --ambiguous=iupac";
$lastzformat="--format=general:name1,zstart1,end1,strand1,name2,zstart2+,end2+,strand2,score,identity,coverage,cigar,cigarx";

$RLIB = "$SCRIPTDIR/rlib";
$ENV{'SCRIPTDIR'} = $SCRIPTDIR;
$ENV{'BCFTOOLS'} = $BCFTOOLS;

# Default parameters (user-modifiable)
$fstem = "contigs";
$kmer = "";
$raw = "";
$wsize = 1000;
$index = "";

# Error messages
$tophelpmsg = "\nUsage:\t snpfilt.pl <command> [options]

Commands:
  pipeline\truns SPAdes assembly, filters assembly and outputs list of SNPs
  snptable\tconstructs table of SNP sites for multiple isolates already analysed by SnpFilt
  plot\t\tproduces PDF showing coverage, MQ and base quality scores near SNP sites

Type 'snpfilt.pl <command>' for more details\n";

%helpmsg = {};
$helpmsg{"pipeline"} = "\nRuns full SnpFilt pipeline generating output files in <outdir> (skipped if existing)
  SPades assembly: contigs.fasta 
  Mpileup file from samtools: contigs.bcf 
  flt_contigs.gz
  masked_contigs.fasta

And, if <refseq.fasta> is supplied
  aln_<fstem>.cigarx
  snpfilt_results_<fstem>.csv
 
FASTQ files are not necessary if pipeline has previously been run and contigs.bam, contigs.bcf or
dat_contigs.gz are already present in <outdir>

Usage:\t snpfilt.pl pipeline <outdir> [reads1.fastq] [reads2.fastq] [refseq.fasta] [options]

Options:
  --fstem -f  filestem of output files [contigs] 
  --kmer  -k  kmer values for SPAdes, comma-deliminated with no spaces [SPAdes default]
  --raw   -r  turn off read-correction before assembly [correction on by default]\n";

$helpmsg{"plot"} = "\nPlots Pdf showing coverage, MQ and base quality scores near variant sites
(i.e. SNP/indels listed in <outdir>/snpfilt_results_<fstem>.csv)
\nUsage:\t snpfilt.pl plot <outdir> <outputpdf> [options]

Options:
  --fstem -f  filestem of output files [contigs] 
  --wsize -w  size of window (on each side) plotted around SNP [1000]
  --index -i  index (1-based) of SNP to be printed, comma-deliminated with no spaces [all]\n";

$helpmsg{"snptable"} = "\nPrints table of SNPs for all isolates in <topdir> (ref base and position in terms of <refseq>)
Note that 'snpfilt.pl pipeline' should have previously been run on all isolates already
  i.e. output file flt_contigs.gz must exist in each directory inside <topdir>
Supplying optional argument <refseq.fasta> will force alignment step to be rerun for all isolates
\nUsage:\t snpfilt.pl snptable <topdir> <outfile.csv> [refseq.fasta] [options]

Options:
  -fstem -f filestem of output files [contigs]\n";


GetOptions("fstem=s" => \$fstem, "kmer=s" => \$kmer, "raw" => \$raw, "wsize=i" => \$wsize,  "index=s" => \$index) or 
    die "Error: invalid options\n";

die "$tophelpmsg\n" if $#ARGV == -1;
$command = $ARGV[0];
if ($#ARGV == 0)
{
    if (exists($helpmsg{$command}))
    {
        die "$helpmsg{$command}\n";
    }
    else
    {
        die "$tophelpmsg"."Error: command \'$ARGV[0]\' unknown\n";
    }
}

if ($command eq "pipeline")
{
    die "$helpmsg{$command}\n" if $#ARGV < 1;
    $outdir = $ARGV[1];
    if ($#ARGV == 2)
    {
        $refseq = $ARGV[2];
    }
    else
    {
        $reads1 = $ARGV[2];
        $reads2 = $ARGV[3];
        $refseq = $ARGV[4] if $#ARGV >= 4;
    }
    $ctgfile = "$outdir/contigs.fasta";
    $fltfile = "$outdir/flt_contigs.gz";
    $datfile = "$outdir/dat_contigs.gz";
    $bcffile = "$outdir/contigs.bcf";
    $bamfile = "$outdir/contigs.bam";
    $bamstem = "$outdir/contigs";                       # suffix is automatically added
    $maskedfasta = "$outdir/masked_contigs.fasta";   
    
    $alnfile="$outdir/aln_$fstem.cigarx";
    $tmpaln="$outdir/$fstem\_tmp.cigarx";
    $tmpsnpfile = "$outdir/tmp_pos.csv";
    $rawsnpfile = "$outdir/rawsnp_$fstem.csv";
    $fltsnpfile = "$outdir/fltsnp_$fstem.csv";
    $resfile = "$outdir/snpfilt_results_$fstem.csv";

    $spadesparam = "";
    $spadesparam = "-k $kmer" if $kmer;
    $spadesparam .= "--only-assembler" if $raw;

    if (!(-e $ctgfile && -s $ctgfile))
    {
        system("$SPADES -o $outdir $spadesparam -1 $reads1 -2 $reads2");
        die "Error: SPAdes assembly exited with error code $?\n" if $?;
    }
    if (!(-e $datfile) || -s $datfile < 100)
    {
        if (!(-e $bcffile) || -s $bcffile < 100)
        {
            if (!(-e $bamfile) || -s $bamfile < 100)
            {
                die "Error: Missing FASTQ files\n" if $#ARGV < 3;
                $samfile = "$outdir/contigs.sam";
                $tmpbamfile = "$outdir/contigs_unsorted.bam";
                
                system ("$BWA index $ctgfile -p $ctgfile.cgi");
                system ("$BWA mem -M $ctgfile $reads1 $reads2 > $samfile");
                system ("$SAMTOOLS faidx $ctgfile");
                system ("$SAMTOOLS view -bt $ctgffile.fai -o $tmpbamfile $samfile");
                system ("$SAMTOOLS sort $tmpbamfile $bamstem");
                system ("rm -r $samfile $tmpbamfile");
            }
            system ("$SAMTOOLS mpileup -gB -Q 0 -t DP,DP4,SP -f $ctgfile $bamfile > $bcffile");
            system ("$BCFTOOLS index -f $bcffile");
        }
        system ("$SCRIPTDIR/format_bcf.pl $bcffile $datfile");      # need bcftools
    }
    system ("R CMD BATCH --no-save --no-restore \"--args RLIB='$RLIB' outfile='$fltfile' datfile='$datfile' wsize=$wsize\" $SCRIPTDIR/filter_contigs.R $outdir/filter_contigs.Rout");
    system ("$SCRIPTDIR/convert_fasta.pl $fltfile > $maskedfasta");

    if ($refseq)
    {
        die "Missing scorefile $lastzscorefile for LASTZ\n" unless (-e $lastzscorefile);
        system ("$LASTZ $refseq\[multiple,nameparse=darkspace\] $ctgfile\[unmask,nameparse=darkspace\] --output=$tmpaln $lastzformat $lastzargs --score=$lastzscorefile");
        system ("$SCRIPTDIR/best_match.pl $tmpaln > $alnfile");
        system ("$SCRIPTDIR/get_seqvar_pos.pl $alnfile $refseq $ctgfile > $rawsnpfile");
        system ("cut -d',' -f 7,8 $rawsnpfile > $tmpsnpfile");
        system ("zgrep -w -f $tmpsnpfile $fltfile > $fltsnpfile");
        system ("$SCRIPTDIR/merge_snptable.pl $rawsnpfile $fltsnpfile > $resfile");
        system ("rm $tmpaln $rawsnpfile $tmpsnpfile $fltsnpfile");
    }
}

if ($command eq "plot")
{
    die "$helpmsg{$command}\n" if $#ARGV != 2;
    $outdir = $ARGV[1];
    $outpdf = $ARGV[2];
    $snpfile = "$outdir/snpfilt_results_$fstem.csv";
    system ("R CMD BATCH --no-save --no-restore \"--args RLIB='$RLIB' outpdf='$outpdf' snpfile='$snpfile' datfile='$outdir/dat_contigs.gz' wsize=$wsize index='$index'\" $SCRIPTDIR/plot_context.R plot_context.Rout");
}

if ($command eq "snptable")
{
    die "$helpmsg{$command}\n" if $#ARGV < 2;
    $topdir = $ARGV[1];
    $snptableout = $ARGV[2];
    if ($#ARGV > 2)
    {
        die "Missing scorefile $lastzscorefile for LASTZ\n" unless (-e $lastzscorefile);
        $refseq = $ARGV[3];
        @dirs = `ls -d $topdir/*/`;

        foreach $dd (@dirs)
        {
            chomp $dd;
            $ctgfile = "$dd/contigs.fasta";
            $fltfile = "$dd/flt_contigs.gz";
            $alnfile="$dd/aln_$fstem.cigarx";
            $tmpaln="$dd/$fstem\_tmp.cigarx";
            $tmpsnpfile = "$dd/tmp_pos.csv";
            $rawsnpfile = "$dd/rawsnp_$fstem.csv";
            $fltsnpfile = "$dd/fltsnp_$fstem.csv";
            $resfile = "$dd/snpfilt_results_$fstem.csv";

            system ("$LASTZ $refseq\[multiple,nameparse=darkspace\] $ctgfile\[unmask,nameparse=darkspace\] --output=$tmpaln $lastzformat $lastzargs --score=$lastzscorefile");
            system ("$SCRIPTDIR/best_match.pl $tmpaln > $alnfile");
            system ("$SCRIPTDIR/get_seqvar_pos.pl $alnfile $refseq $ctgfile > $rawsnpfile");    # need bcftools
            system ("cut -d',' -f 7,8 $rawsnpfile > $tmpsnpfile");
            system ("zgrep -w -f $tmpsnpfile $fltfile > $fltsnpfile");
            system ("$SCRIPTDIR/merge_snptable.pl $rawsnpfile $fltsnpfile > $resfile");
            system ("rm $tmpaln $rawsnpfile $tmpsnpfile $fltsnpfile");
        }
    }
    system ("$SCRIPTDIR/get_SNPtable.pl $fstem $topdir > $snptableout");  # need SCRIPTDIR
} 
