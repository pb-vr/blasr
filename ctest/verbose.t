Set up
  $ . $TESTDIR/setup.sh

Test alignment score
  $ $EXEC $DATDIR/lambda_bax.fofn  $DATDIR/lambda_ref.fasta -holeNumbers 1-200 -V 3 > $TMP1
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
