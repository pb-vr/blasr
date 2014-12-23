Set up
  $ . $TESTDIR/setup.sh

bug_25766, added an option -minRawSubreadScore
  $ BASFILE=$DATDIR/lambda_bax.fofn
  $ REF=$DATDIR/lambda_ref.fasta
  $ $EXEC $BASFILE $REF -out $TMP1 -minRawSubreadScore 700 -nproc 18
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
  $ sort $TMP1 > $TMP2
  $ diff $TMP2 $STDDIR/bug_25766.m4
