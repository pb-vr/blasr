Set up directories
  $ CURDIR=$TESTDIR
  $ REMOTEDIR=/mnt/secondary-siv/testdata/BlasrTestData/ctest
  $ DATDIR=$REMOTEDIR/data
  $ OUTDIR=$CURDIR/out
  $ STDDIR=$REMOTEDIR/stdout

Set up the executable: sawriter.
  $ EXEC=$TESTDIR/../sawriter

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ $EXEC $OUTDIR/ecoli_larsson.sa $DATDIR/ecoli_reference.fasta -blt 11 -larsson
  $ echo $?
  0

  $ $EXEC $OUTDIR/ecoli_welter.sa $DATDIR/ecoli_reference.fasta -blt 11 -welter 2>$OUTDIR/sawriter.log
  $ echo $?
  0

  $ md5sum $OUTDIR/ecoli_larsson.sa |cut -f 1 -d ' ' 
  e23b6afe6ddd74b2656e36bf93f6840c

  $ md5sum $OUTDIR/ecoli_welter.sa |cut -f 1 -d ' '
  e23b6afe6ddd74b2656e36bf93f6840c
