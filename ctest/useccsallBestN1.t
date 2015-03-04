Set up
  $ . $TESTDIR/setup.sh

Test -useccsall with bestn = 1
  $ $EXEC $DATDIR/ccstest.fofn $DATDIR/ccstest_ref.fasta -bestn 1 -useccsall -sam -out $OUTDIR/useccsall.sam -holeNumbers 76772
  [INFO]* (glob)
  [INFO]* (glob)
  $ sed -n '9,$ p' $OUTDIR/useccsall.sam > $TMP1
  $ sed -n '9,$ p' $STDDIR/useccsall_2015_03_01.sam > $TMP2
# changelist 146599, 2015_03_01, reads from the same SMRTCell should share unique read group ID.
  $ diff $TMP1 $TMP2
  $ rm $TMP1 $TMP2
# useccsall_2014_08_21.sam --> change 138516: added YS, YE, ZM tags
# useccsall_2014_08_28.sam --> change 139176: Update SAM MD5 
