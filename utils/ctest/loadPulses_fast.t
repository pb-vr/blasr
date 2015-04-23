Set up directories
  $ . $TESTDIR/setup.sh 

Set up the executable: loadPulses.
  $ EXEC=$TESTDIR/../loadPulses

#Test loadPulses: input is a bas.h5 file
#Test -byread and -bymetric
  $ BAS_IN_1=$DATDIR/lambda_lp.fofn
  $ CMP_IN_1=$DATDIR/lambda_lp.cmp.h5
  $ CMP_STDOUT_1=$STDDIR/lambda_lp.cmp.h5
  $ CMP_OUT_byread_1=$OUTDIR/lambda.byread.cmp.h5
  $ CMP_OUT_bymetric_1=$OUTDIR/lambda.bymetric.cmp.h5
  $ METRICS=QualityValue,MergeQV,InsertionQV,DeletionQV,DeletionTag,PulseWidth,SubstitutionQV,SubstitutionTag

  $ rm -f $CMP_OUT_byread_1
  $ cp $CMP_IN_1 $CMP_OUT_byread_1
  $ $EXEC $BAS_IN_1 $CMP_OUT_byread_1 -metrics $METRICS -byread > $OUTDIR/tmp.log
  [INFO] * [loadPulses] started. (glob)
  [INFO] * [loadPulses] ended. (glob)

  $ h5diff -c $CMP_OUT_byread_1 $CMP_STDOUT_1
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

  $ rm -f $CMP_OUT_bymetric_1
  $ cp $CMP_IN_1 $CMP_OUT_bymetric_1
  $ $EXEC $BAS_IN_1 $CMP_OUT_bymetric_1 -metrics $METRICS -bymetric > $OUTDIR/tmp.log 
  [INFO] * [loadPulses] started. (glob)
  [INFO] * [loadPulses] ended. (glob)
  $ h5diff -c $CMP_OUT_bymetric_1 $CMP_STDOUT_1 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]


