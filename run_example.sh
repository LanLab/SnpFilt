# This script runs snpfilt on the sample data sets in the examples directory
# Comments are provided to show different forms of usage
# The script can be run in batch mode but note that some commands can take a long time

cd examples

# Modify this if snpfilt was not installed in the default directory 
SFPATH=../snpfilt

# Analysis takes about an hour including assembly
gunzip L1701_S15_L001_R1_001.fastq.gz
gunzip L1701_S15_L001_R2_001.fastq.gz
$SFPATH/snpfilt.pl pipeline L1701 L1701_S15_L001_R1_001.fastq L1701_S15_L001_R2_001.fastq lt2_refchrom.fasta

# Re-analysis is faster when assembly and intermediate files do not need to regenerated
# The '-f' flag is used to specify output file names otherwise previous files are overwritten
# Note that for re-analysis, fastq files are no longer required
# If space is constrained, you can delete fastq files, contigs.bam and even contig.bcf, although
# the third file is considerably smaller than the others and may be useful other applications)
$SFPATH/snpfilt.pl pipeline L1701 dt104_refchrom.fasta -f rDT104

# Assembly and filtering can be performed without any reference sequence
gunzip L1688_S14_L001_R1_001.fastq.gz
gunzip L1688_S14_L001_R2_001.fastq.gz
$SFPATH/snpfilt.pl pipeline L1688 L1688_S14_L001_R1_001.fastq L1688_S14_L001_R2_001.fastq

# Visually examine contig region around SNPs
# (But not for L1688 and there are no SNPs without a reference provided)
$SFPATH/snpfilt.pl plot L1701 L1701_snpAll.pdf

# Examine only SNPS 10, 20 and 30 with a larger neighbouring region
$SFPATH/snpfilt.pl plot L1701 L1701_i10i20i30_w2000.pdf -i 10,20,30 -w 2000

# Output table of SNPs for L1701 and L1688
$SFPATH/snpfilt.pl pipeline L1688 lt2_refchrom.fasta
$SFPATH/snpfilt.pl snptable . lt2_snptable.csv

# The following command will not run as L1688 has not been aligned to dt104 
# SFPATH/snpfilt.pl snptable . dt104_snptable.csv -f rDT104
# But snpfilt will automatically rerun alignment if reference sequence is provided
$SFPATH/snpfilt.pl snptable . dt104_snptable.csv dt104_refchrom.fasta -f rDT104

