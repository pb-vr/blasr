Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: sa2bwt.
  $ EXEC=$TESTDIR/../sa2bwt

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ FA=$DATDIR/ecoli_reference.fasta
  $ SA=$DATDIR/ecoli_reference.sa
  $ BWT=$OUTDIR/ecoli_reference.bwt
  $ $EXEC $FA $SA $BWT 
  $ echo $?
  0
