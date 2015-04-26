# Set up directories
CURDIR=$TESTDIR
REMOTEDIR=/mnt/secondary-siv/testdata/BlasrTestData/ctest
DATDIR=$REMOTEDIR/data
OUTDIR=$CURDIR/out
STDDIR=$REMOTEDIR/stdout

# Set up the executable: blasr.
EXEC=$TESTDIR/../blasr

# Define tmporary files
TMP1=$OUTDIR/$$.tmp.out
TMP2=$OUTDIR/$$.tmp.stdout

# Make OUTDIR
mkdir -p $OUTDIR

#FIXME: make samtools independent of absolute build path.
SAMTOOLS=/mnt/secondary/Smrtpipe/builds/Internal_Mainline_Nightly_LastSuccessfulBuild/analysis/bin/samtools

#Update date
UPDATEDATE=20150425 # changelist, update CIGAR string, replace M with X=
