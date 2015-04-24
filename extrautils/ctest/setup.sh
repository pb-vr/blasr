# Set up directories
CURDIR=$TESTDIR
REMOTEDIR=/mnt/secondary-siv/testdata/BlasrTestData/ctest
DATDIR=$REMOTEDIR/data
OUTDIR=$CURDIR/out
STDDIR=$REMOTEDIR/stdout

# Define tmporary files
TMP1=$OUTDIR/$$.tmp.out
TMP2=$OUTDIR/$$.tmp.stdout

# Make OUTDIR
mkdir -p $OUTDIR

