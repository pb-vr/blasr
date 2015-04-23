Set up directories
  $ . $TESTDIR/setup.sh 

Set up the executable: loadPulses.
  $ EXEC=$TESTDIR/../loadPulses

#Test loadPulses: input is a pls.h5 file
#Test -byread and -bymetric
  $ PLS_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN_2=$DATDIR/ecoli_lp_tiny.cmp.h5
  $ CMP_STDOUT_2=$STDDIR/ecoli_lp_tiny.cmp.h5
  $ CMP_OUT_byread_2=$OUTDIR/ecoli_lp_tiny.byread.cmp.h5
  $ CMP_OUT_bymetric_2=$OUTDIR/ecoli_lp_tiny.bymetric.cmp.h5
  $ METRICS=StartFrame,PulseWidth,WidthInFrames,pkmid,IPD,Light

  $ rm -f  $CMP_OUT_byread_2
  $ cp $CMP_IN_2 $CMP_OUT_byread_2
  $ $EXEC $PLS_IN $CMP_OUT_byread_2 -metrics $METRICS -byread
  [INFO] * [loadPulses] started. (glob)
  loading 2 alignments for movie 1
  loading 2 alignments for movie 2
  [INFO] * [loadPulses] ended. (glob)

  $ h5diff -c $CMP_OUT_byread_2 $CMP_STDOUT_2
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

  $ rm -f $CMP_OUT_bymetric_2
  $ cp $CMP_IN_2 $CMP_OUT_bymetric_2
  $ $EXEC $PLS_IN $CMP_OUT_bymetric_2 -metrics $METRICS -bymetric 
  [INFO] * [loadPulses] started. (glob)
  loading 2 alignments for movie 1
  loading 2 alignments for movie 2
  [INFO] * [loadPulses] ended. (glob)

  $ h5diff -c $CMP_OUT_bymetric_2 $CMP_STDOUT_2 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

#Test loadPulses for deep sorted cmp.h5
  $ FOFN_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN_SORTED=$DATDIR/ecoli_lp_tiny_sorted.cmp.h5
  $ CMP_STDOUT_SORTED=$STDDIR/ecoli_lp_tiny_sorted.cmp.h5
  $ CMP_OUT_SORTED_bymetric=$OUTDIR/ecoli_lp_tiny_sorted_bymetric.cmp.h5
  $ CMP_OUT_SORTED_byread=$OUTDIR/ecoli_lp_tiny_sorted_byread.cmp.h5
  $ METRICS=StartFrame,PulseWidth,WidthInFrames,pkmid,IPD,Light,DeletionQV,InsertionQV,SubstitutionQV,MergeQV,QualityValue,DeletionTag,SubstitutionTag,ClassifierQV,PreBaseFrames,PulseIndex

  $ rm -f $CMP_OUT_SORTED_bymetric
  $ cp $CMP_IN_SORTED $CMP_OUT_SORTED_bymetric
  $ $EXEC $FOFN_IN $CMP_OUT_SORTED_bymetric -bymetric -metrics $METRICS > $OUTDIR/tmp.log
  [INFO] * [loadPulses] started. (glob)
  [INFO] * [loadPulses] ended. (glob)

  $ h5diff -c $CMP_OUT_SORTED_bymetric $CMP_STDOUT_SORTED 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]


  $ rm -f  $CMP_OUT_SORTED_byread
  $ cp $CMP_IN_SORTED $CMP_OUT_SORTED_byread
  $ $EXEC $FOFN_IN $CMP_OUT_SORTED_byread -byread -metrics $METRICS > $OUTDIR/tmp.log
  [INFO] * [loadPulses] started. (glob)
  [INFO] * [loadPulses] ended. (glob)

  $ h5diff -c $CMP_OUT_SORTED_bymetric $CMP_STDOUT_SORTED 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

#Test loadPulses for a zero-alignment cmp.h5 file.
  $ FOFN_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN_NOALN=$DATDIR/noaln_lp.cmp.h5
  $ $EXEC $FOFN_IN $CMP_IN_NOALN -byread -metrics $METRICS
  [INFO] * [loadPulses] started. (glob)
  WARNING, there is no alignment in the cmp file.
  [INFO] * [loadPulses] ended. (glob)

#Test loadPulses -byMetric with a 'large' bas.h5 file of which the dataset size is greater than maxElements. 
  $ FOFN_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN=$DATDIR/ecoli_lp_tiny_sorted.cmp.h5
  $ CMP_OUT=$OUTDIR/ecoli_lp_maxEle.cmp.h5
  $ CMP_STDOUT=$STDDIR/ecoli_lp_maxEle.cmp.h5
  $ METRICS=QualityValue,MergeQV,InsertionQV,DeletionQV,DeletionTag,PulseWidth,SubstitutionQV,SubstitutionTag
  $ MAX_ELEMENTS=140000000

  $ rm -f $CMP_OUT
  $ cp $CMP_IN $CMP_OUT
  $ $EXEC $FOFN_IN $CMP_OUT -bymetric -metrics $METRICS -maxElements $MAX_ELEMENTS
  [INFO] * [loadPulses] started. (glob)
  Either the number of elements exceeds maxElement (140000000). Or the estimated memory 
  consumption exceeds maxMemory (4 GB).
  Loading pulses from .+ by read. (re)
  loading 2 alignments for movie 1
  loading 2 alignments for movie 2
  [INFO] * [loadPulses] ended. (glob)

  $ h5diff -c $CMP_OUT $CMP_STDOUT
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  dataset: </FileLog/Version> and </FileLog/Version>
  \d+ differences found (re)
  [1]

#Test loadPulses -byMetric on a multi-streaming job.
  $ FOFN_IN=$DATDIR/lambda_bax.fofn
  $ CMP_IN=$DATDIR/lambda_bax.cmp.h5
  $ CMP_OUT=$OUTDIR/lambda_bax.cmp.h5
  $ CMP_STDOUT=$STDDIR/lambda_bax.cmp.h5
  $ METRICS=QualityValue,MergeQV,InsertionQV,DeletionQV,DeletionTag,PulseWidth,SubstitutionQV,SubstitutionTag

  $ rm -f $CMP_OUT
  $ cp $CMP_IN $CMP_OUT
  $ $EXEC $FOFN_IN $CMP_OUT -bymetric -metrics $METRICS 
  [INFO] * [loadPulses] started. (glob)
  WARNING: There is insufficient data to compute metric: MergeQV in the file .+ It will be ignored. (re)
  loading 2 alignments for movie 1
  loading 2 alignments for movie 1
  [INFO] * [loadPulses] ended. (glob)

  $ h5diff -c $CMP_OUT $CMP_STDOUT
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  dataset: </FileLog/Version> and </FileLog/Version>
  \d+ differences found (re)
  [1]


#Test loadPulses -bymetric for a ccs cmp.h5 file generated from multiple movies.
  $ FOFN_IN=$DATDIR/ccs_lp.fofn
  $ CMP_IN=$DATDIR/ccs_lp.cmp.h5
  $ CMP_OUT=$OUTDIR/ccs_lp.cmp.h5

# The original pls.h5 files disappeared, to use another dataset instead.
  $ rm -f CMP_OUT
  $ cp $CMP_IN $CMP_OUT
  $ $EXEC $FOFN_IN $CMP_OUT -bymetric -metrics QualityValue
  [INFO] * [loadPulses] started. (glob)
  loading 100 alignments for movie 1
  loading 45 alignments for movie 2
  [INFO] * [loadPulses] ended. (glob)

  $ h5ls -r $CMP_OUT | grep "AlnInfo"
  /AlnInfo                 Group
  /AlnInfo/AlnIndex        Dataset {145/Inf, 22}
  /AlnInfo/NumPasses       Dataset {145/Inf}


#Test loadPulses *.fofn cmp.h5 where *.fofn can either contain ccs.h5 or bas.h5
# and the cmp.h5 s readType is CCS
  $ CCS_FOFN=$DATDIR/test_ccs.fofn
  $ BAS_FOFN=$DATDIR/test_bas.fofn

  $ CMP_IN=$DATDIR/test_ccs_bas.cmp.h5
  $ CCS_OUT=$OUTDIR/test_ccs_bas_ccs.cmp.h5
  $ BAS_OUT=$OUTDIR/test_ccs_bas_bas.cmp.h5

  $ cp $CMP_IN $CCS_OUT 
  $ cp $CMP_IN $BAS_OUT

  $ $EXEC $CCS_FOFN $CCS_OUT -metrics QualityValue,DeletionQV,DeletionTag,InsertionQV,SubstitutionQV
  [INFO] * [loadPulses] started. (glob)
  loading 11 alignments for movie 1
  [INFO] * [loadPulses] ended. (glob)

  $ sleep 1

  $ $EXEC $BAS_FOFN $BAS_OUT -metrics QualityValue,DeletionQV,DeletionTag,InsertionQV,SubstitutionQV
  [INFO] * [loadPulses] started. (glob)
  loading 11 alignments for movie 1
  [INFO] * [loadPulses] ended. (glob)
 
  $ h5diff $CCS_OUT $BAS_OUT
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine> (glob)
  * differences found (glob)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  * differences found (glob)
  [1]

