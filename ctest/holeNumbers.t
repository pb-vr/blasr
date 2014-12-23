Set up
  $ . $TESTDIR/setup.sh

Test -holeNumbers
  $ rm -f $OUTDIR/holeNumbers.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out $OUTDIR/holeNumbers.m4 -holeNumbers 14798,55000-55100 -nproc 8
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/holeNumbers.m4 > $TMP1
  $ sort $STDDIR/holeNumbers_2014_05_29.m4 > $TMP2
  $ diff $TMP1 $TMP2
  $ rm $TMP1 $TMP2
