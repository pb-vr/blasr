Set up
  $ . $TESTDIR/setup.sh

Test blasr with *.fofn input
#  $ rm -rf $OUTDIR/lambda_bax.m4 
#  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 --out lambda_bax_tmp.m4 --nproc 15 --minMatch 14
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ sort lambda_bax_tmp.m4 > $OUTDIR/lambda_bax.m4
#  $ diff $OUTDIR/lambda_bax.m4 $STDDIR/lambda_bax.m4
# This test takes a long time, use a subset instad. 

  $ rm -rf $OUTDIR/lambda_bax_subset.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 --out $OUTDIR/lambda_bax_tmp_subset.m4 --nproc 15 --minMatch 14 --holeNumbers 1--1000 --sa $DATDIR/lambda_ref.sa
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/lambda_bax_tmp_subset.m4 > $OUTDIR/lambda_bax_subset.m4
  $ diff $OUTDIR/lambda_bax_subset.m4 $STDDIR/lambda_bax_subset.m4
