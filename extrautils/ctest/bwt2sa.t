Set up 
  $ . $TESTDIR/setup.sh

Noticed that pipeline fa->sa->bwt works OK, however fa->sa->bwt->sa 
does not generate identical suffix array.
Set up the executable: bwt2sa.
  $ EXEC=$TESTDIR/../bwt2sa

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ SA=$OUTDIR/ecoli_reference.bwt2sa.sa
  $ BWT=$DATDIR/ecoli_reference.bwt
  $ $EXEC $BWT $SA
  $ echo $?
  0

