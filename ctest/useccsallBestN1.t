Set up
  $ . $TESTDIR/setup.sh

Test -useccsall with bestn = 1
  $ $EXEC $DATDIR/ccstest.fofn $DATDIR/ccstest_ref.fasta -bestn 1 -useccsall -sam -out $OUTDIR/useccsall.sam -holeNumbers 76772
  [INFO]* (glob)
  [INFO]* (glob)
  $ sed -n '9,$ p' $OUTDIR/useccsall.sam > $TMP1
  $ sed -n '9,$ p' $STDDIR/$UPDATEDATE/useccsall.sam > $TMP2
  $ diff $TMP1 $TMP2
  $ rm $TMP1 $TMP2
