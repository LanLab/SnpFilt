#!/bin/bash
CURRDIR=$(pwd)

################
# User-modifiable variables

# Directory where the program is installed
#SCRIPTDIR="$CURRDIR/snpfilt"
SCRIPTDIR=/srv/scratch/lanlab/SnpFilt

# Dependent programs
# By default, the script will look for versions on the path
# Paths can be manually changed by user before running install.sh
# If they are not found (or invalid path supplied) then a new version  be installed inside SCRIPTDIR 
# (source code inside third_party)

SAMTOOLS="samtools"     # look for system version (will install into snpfilt directory if not found)
BCFTOOLS="bcftools"
SPADES="spades.py"
LASTZ="lastz"
BWA="bwa"

#################
# Starting script

for rqcmd in perl R zgrep gunzip cut
do
    echo -n "Looking for $rqcmd"
    res=$(type -P $rqcmd)
    if [ -n $res ]
    then
        echo " ... $res [OK]"
    else
        echo "Didn't find $rqcmd. Please install and put in path before continuing"
        exit 1
    fi
done

for module in Getopt::Long
do
    echo -n "Looking for $module"
    res=$(perldoc -lm "Getopt::Long")
    if [ -n $res ]
    then
        echo " ... $res [OK]"
    else
        echo "    Module $module not found. Please install before continuing"
        exit 1
    fi
done

if [ ! -d $SCRIPTDIR ]
then
    mkdir -p $SCRIPTDIR
fi

echo -n "Checking LASTZ ..."
res=$(type -P $LASTZ)
if [ -z $res ]
then 
    echo -n "installing into $SCRIPTDIR"
    cd third_party
    tar -xvzf lastz-1.02.00.tar.gz
    cd lastz-distrib-1.02.00/src
    export LASTZ_INSTALL=$SCRIPTDIR
    make
    make install
    if [ $? -eq 0 ]
    then
        LASTZ="${SCRIPTDIR}/lastz"
        echo " [OK]"
        cd $CURRDIR
        rm -r third_party/lastz-distrib-1.02.00
    else
        echo "Could not install LASTZ"
        exit 1
    fi
    unset LASTZ_INSTALL
else
    echo " $LASTZ [OK]"
fi

echo -n "Checking samtools ..."
res=$(type -P $SAMTOOLS)
if [ -z $res ]
then 
    echo "installing into $SCRIPTDIR"
    cd third_party
    tar xvjf samtools-1.2.tar.bz2
    cd samtools-1.2
    make
    make prefix=$SCRIPTDIR install
    if [ $? -eq 0 ]
    then
        SAMTOOLS="${SCRIPTDIR}/bin/samtools"
        echo " [OK]"
        cd $CURRDIR
        rm -r third_party/samtools-1.2
    else
        echo "Could not install samtools"
        exit 1
    fi
else
    SAMTOOLS=$res
    echo " $SAMTOOLS [OK]"
fi

echo -n "Checking bcftools ..."
res=$(type -P $BCFTOOLS)
if [ -z $res ]
then 
    echo "installing into $SCRIPTDIR"
    cd third_party
    tar xvjf bcftools-1.2.tar.bz2
    cd bcftools-1.2
    make
    make prefix=$SCRIPTDIR install
    if [ $? -eq 0 ]
    then
        BCFTOOLS="${SCRIPTDIR}/bin/bcftools"
        echo " [OK]"
        cd $CURRDIR
        rm -r third_party/bcftools-1.2
    else
        echo "Could not install samtools"
        exit 1
    fi
else
    BCFTOOLS=$res
    echo " $BCFTOOLS [OK]"
fi

echo -n "Checking spades ..."
res=$(type -P $SPADES)
if [ -z $res ]
then 
    echo "installing into $SCRIPTDIR"
    cd third_party
    tar -xvzf SPAdes-3.6.1-Linux.tar.gz
    mv SPAdes-3.6.1-Linux $SCRIPTDIR
    SPADES="${SCRIPTDIR}/SPAdes-3.6.1-Linux/bin/spades.py"
    cd $CURRDIR
else
    SPADES=$res
    echo " $SPADES [OK]"
fi

echo -n "Checking BWA ..."
res=$(type -P $BWA)
if [ -z $res ]
then 
    echo "installing into $SCRIPTDIR"
    cd third_party
    tar -xvjf bwa-0.7.12.tar.bz2
    mv bwa-0.7.12 $SCRIPTDIR
    cd ${SCRIPTDIR}/bwa-0.7.12
    make
    BWA="${SCRIPTDIR}/bwa-0.7.12/bwa"
    cd $CURRDIR
else
    BWA=$res
    echo " $BWA [OK]"
fi

# Install required R-package inside SnpFilt directory
# We don't use default user libraries to avoid problems
# when installing for in multi-user environments
mkdir $SCRIPTDIR/rlib
R CMD BATCH --vanilla "--args RLIB='$SCRIPTDIR/rlib'" install_rpckg.R
if [ $? -ne 0 ]
then
    echo "error: could not install R package - check install_rpckg.Rout"
    exit 1
fi

cd src
cp best_match.pl convert_fasta.pl filter_contigs.R find_ctgpos.pl format_bcf.pl get_seqvar_pos.pl get_SNPtable.pl *.scores merge_snptable.pl plot_context.R $SCRIPTDIR
outscript="$SCRIPTDIR/snpfilt.pl"
echo "#!/bin/env perl" > $outscript
echo "" >> $outscript
echo "\$SCRIPTDIR = \"${SCRIPTDIR}\";"  >> $outscript
echo "\$SAMTOOLS = \"${SAMTOOLS}\";" >> $outscript
echo "\$BCFTOOLS = \"${BCFTOOLS}\";" >> $outscript
echo "\$LASTZ = \"${LASTZ}\";" >> $outscript
echo "\$SPADES = \"${SPADES}\";" >> $outscript
echo "\$BWA = \"${BWA}\";" >> $outscript
cat snpfilt_stem.pl >> $outscript

chmod a+x $SCRIPTDIR
chmod a+x $SCRIPTDIR/snpfilt.pl

echo ""
echo "Install complete!"
echo "To place snpfilt into PATH, append $SCRIPTDIR into path or symlink to $SCRIPTDIR/snpfilt.pl"
