Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: sals.
  $ EXEC=$TESTDIR/../sals

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ $EXEC $DATDIR/ecoli_reference.sa 
   * has a suffix array.
   * has a lookup table for word size. 8
  $ echo $?
  0
