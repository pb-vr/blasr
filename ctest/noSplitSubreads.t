Set up
  $ . $TESTDIR/setup.sh

Test blasr with --noSplitSubreads 
#  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads.m4 
#  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta --noSplitSubreads -m 4 --out lambda_bax_noSplitSubreads_tmp.m4 --nproc 15
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ sort lambda_bax_noSplitSubreads_tmp.m4 > $OUTDIR/lambda_bax_noSplitSubreads.m4
#  $ diff $OUTDIR/lambda_bax_noSplitSubreads.m4 $STDDIR/lambda_bax_noSplitSubreads.m4
# This test takes a long time, use a subset instad. 

  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads_subset.m4 
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta --noSplitSubreads -m 4 --out $OUTDIR/lambda_bax_noSplitSubreads_tmp_subset.m4 --nproc 15 --holeNumbers 1--1000 --sa $DATDIR/lambda_ref.sa
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/lambda_bax_noSplitSubreads_tmp_subset.m4 > $OUTDIR/lambda_bax_noSplitSubreads_subset.m4
  $ diff $OUTDIR/lambda_bax_noSplitSubreads_subset.m4 $STDDIR/lambda_bax_noSplitSubreads_subset.m4
