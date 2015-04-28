Set up directories
  $ CURDIR=$TESTDIR
  $ REMOTEDIR=/mnt/secondary-siv/testdata/BlasrTestData/ctest
  $ DATDIR=$REMOTEDIR/data
  $ OUTDIR=$CURDIR/out
  $ STDDIR=$REMOTEDIR/stdout

Set up the executable: sdpMatcher.
  $ EXEC=$TESTDIR/../sdpMatcher

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ FA=$DATDIR/ecoli_subset.fasta
  $ $EXEC $FA $FA 10 -local  > $OUTDIR/sdpMatcher.out
  $ echo $?
  0
  $ diff $OUTDIR/sdpMatcher.out $STDDIR/sdpMatcher.stdout


