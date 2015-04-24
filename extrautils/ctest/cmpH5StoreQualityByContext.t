Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: cmpH5StoreQualityByContext.
  $ EXEC=$TESTDIR/../cmpH5StoreQualityByContext

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ $EXEC $DATDIR/ecoli_out.cmp.h5 $OUTDIR/ecoli_out.qbc -contextLength 8 -onlyMaxLength -minSamples 600 -maxSamples 1500 > $TMP1
  $ echo $?
  0
