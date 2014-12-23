Set up
  $ . $TESTDIR/setup.sh

Test blasr with -m 0 ~ 5 
  $ rm -rf $OUTDIR/read.m0
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 0 -out $OUTDIR/read.m0
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/read.m0 $STDDIR/read.m0

  $ rm -rf $OUTDIR/read.m1
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 1 -out $OUTDIR/read.m1
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/read.m1 $STDDIR/read_2014_05_29.m1

  $ rm -rf $OUTDIR/read.m2
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 2 -out $OUTDIR/read.m2
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/read.m2 $STDDIR/read.m2

  $ rm -rf $OUTDIR/read.m3
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 3 -out $OUTDIR/read.m3
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/read.m3 $STDDIR/read.m3

  $ rm -rf $OUTDIR/read.m4
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 4 -out $OUTDIR/read.m4
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/read.m4 $STDDIR/read.m4
