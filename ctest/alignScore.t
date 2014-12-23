Set up
  $ . $TESTDIR/setup.sh

Test alignment score
  $ rm -rf $OUTDIR/testscore.m0
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -minReadLength 1 -m 0 -out $OUTDIR/testscore.m0
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/testscore.m0 $STDDIR/testscore.m0
