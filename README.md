SnpFilt is a package for identifying sequence variants in high-throughput
sequencing reads using reference-free assemblies

Installation
------------
SnpFilt has only been tested on a Linux system. Specifically it requires the
following programs to be available
* Perl
* R
* gunzip
* zgrep
* cut

It also requires a number of third-party packages
* samtools
* bcftools
* BWA
* SPAdes
* LASTZ
* caTools (R-package)

By default, install.sh will check if any of these can be found in the path
If not it will install missing pacakges from the third_party directory
The user can direct the install script to other locations by modifying
lines 16-20 in install.sh. Set variable to empty string to force re-installation

By default, the SnpFilt package will be installed in a folder called
snpfilt which will be created inside the current directory,
but the user can change this by modifying line 8 in install.sh

Once file paths are set correctly, SnpFilt can be installed typing

    > ./install.sh

To put in snpfilt into the path you can add the snpfilt directory to your $PATH varable
or put a symlink to snpfilt/snpfilt.pl in a directory in your path.


Running SnpFilt
---------------
To find sequence variants to a given set of reads type

    > <snpfiltpath>/snpfilt.pl pipeline <outdir> <reads1.fastq> <reads2.fastq> <refseq.fasta>

or 

    > perl <snpfiltpath>/snpfilt.pl pipeline <outdir> <reads1.fastq> <reads2.fastq> <refseq.fasta>

For further instructions regarding usage simply type

    > snpfilt.pl 

Test data is provided in the examples directory and detailed comments 
regarding usage is available in the run_example.sh script.


Output Files
------------
SnpFilt creates a directory for each isolate and generates a number of
output files in the directory

* contigs.fasta: this the raw assembly generated by SPAdes, which is likely to contain sequence errors
* masked_contigs.fasta: this is assembly with problematic sites masked to 'N'
* flt_contigs.gz: compressed CSV file listing the base and filter code for each base in the assembly. Only sites with flt=0 are reliable. See snpfilt_results_contigs.csv for an explanation on the other filter codes.
* contigs.bcf: mpileup file generated by mapping the read back to contigs.fasta

The files listed above are always created and use only information based on
the reads themself. When a reference genome is provided the isolate will be
compared with it generate the additional files
* snpfilt_results_fstem.csv: CSV file listing all potential SNPs and the filter code describing their reliability
* aln_fstem.cigarx: tabular file giving positions where assembly aligned to the reference genome



