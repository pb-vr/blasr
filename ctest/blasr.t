Set up directories
  $ CURDIR=$TESTDIR
  $ REMOTEDIR=/mnt/secondary-siv/secondarytest/testdata/BlasrTestData/ctest
  $ DATDIR=$REMOTEDIR/data
  $ OUTDIR=$CURDIR/out
  $ STDDIR=$REMOTEDIR/stdout

Set up the executable: blasr.
  $ BIN=$TESTDIR/../alignment/bin
  $ EXEC=$BIN/blasr

Test blasr on ecoli.
Test blasr with -sam
#See $STDOUT/ecoli.sam for 1.4 output.
  $ rm -rf $OUTDIR/ecoli.sam
  $ $EXEC $DATDIR/ecoli.fasta $DATDIR/ecoli_reference.fasta -sam -out $OUTDIR/ecoli.sam -nproc 15
  $ tail -n+5 $OUTDIR/ecoli.sam | sort | cut -f 1-11| md5sum
  d2e2b6cfe710b7b6a065e73c3244b0b9  -

#  $ tail -n+5 $STDDIR/ecoli.sam | sort | cut -f 1-11| md5sum
#  d2e2b6cfe710b7b6a065e73c3244b0b9  -

Test blasr with -m 0 ~ 5 
  $ rm -rf $OUTDIR/read.m0
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 0 -out $OUTDIR/read.m0
  $ diff $OUTDIR/read.m0 $STDDIR/read.m0

  $ rm -rf $OUTDIR/read.m1
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 1 -out $OUTDIR/read.m1
  $ diff $OUTDIR/read.m1 $STDDIR/read.m1

  $ rm -rf $OUTDIR/read.m2
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 2 -out $OUTDIR/read.m2
  $ diff $OUTDIR/read.m2 $STDDIR/read.m2

  $ rm -rf $OUTDIR/read.m3
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 3 -out $OUTDIR/read.m3
  $ diff $OUTDIR/read.m3 $STDDIR/read.m3

  $ rm -rf $OUTDIR/read.m4
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 4 -out $OUTDIR/read.m4
  $ diff $OUTDIR/read.m4 $STDDIR/read.m4

Test blasr with *.fofn input
  $ rm -rf $OUTDIR/lambda_bax.m4 
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out lambda_bax_tmp.m4 -nproc 15
  $ sort lambda_bax_tmp.m4 > $OUTDIR/lambda_bax.m4
  $ diff $OUTDIR/lambda_bax.m4 $STDDIR/lambda_bax.m4

Test blasr with -noSplitSubreads 
  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads.m4 $OUTDIR/lambda_bax_noSplitSubreads.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -noSplitSubreads -m 4 -out lambda_bax_noSplitSubreads_tmp.m4 -nproc 15
  $ sort lambda_bax_noSplitSubreads_tmp.m4 > $OUTDIR/lambda_bax_noSplitSubreads.m4
  $ diff $OUTDIR/lambda_bax_noSplitSubreads.m4 $STDDIR/lambda_bax_noSplitSubreads.m4

Test alignment score
  $ rm -rf $OUTDIR/testscore.m0
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -minReadLength 1 -m 0 -out $OUTDIR/testscore.m0
  $ diff $OUTDIR/testscore.m0 $STDDIR/testscore.m0

Test affineAlign
  $ rm -rf $OUTDIR/affineAlign.m0
  $ $EXEC $DATDIR/affineAlign.fofn $DATDIR/substr_with_ins.fasta -m 0 -out $OUTDIR/affineAlign.m0  -affineAlign  -readIndex 493 -insertion 100 -deletion 100
  $ diff $OUTDIR/affineAlign.m0 $STDDIR/affineAlign.m0

  $ rm -rf $OUTDIR/ecoli_affine.m0
  $ $EXEC $DATDIR/ecoli_affine.fasta $DATDIR/ecoli_reference.fasta -m 0 -out $OUTDIR/ecoli_affine.m0 -affineAlign -insertion 100 -deletion 100
  $ diff $OUTDIR/ecoli_affine.m0 $STDDIR/ecoli_affine.m0


Test -holeNumbers
  $ rm -f $OUTDIR/holeNumbers.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out $OUTDIR/holeNumbers.m4 -holeNumbers 14798,55000-55100 -nproc 8
  $ sort $OUTDIR/holeNumbers.m4 | md5sum
  21fd37b14b85ef7dda332ea10edc524a  -

Test Sam out nm tag
  $ rm -rf $OUTDIR/read.sam
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -sam -out $OUTDIR/read.sam
  $ tail -n+5 $OUTDIR/read.sam |cut -f 18 
  NM:i:2
  NM:i:3
  NM:i:2
  NM:i:4

Test -useccsall with bestn = 1
  $ $EXEC $DATDIR/ccstest.fofn $DATDIR/ccstest_ref.fasta -bestn 1 -useccsall -sam -out $OUTDIR/useccsall.sam -holeNumbers 76772
  $ tail -n+9 $OUTDIR/useccsall.sam | md5sum 
  45fa5f07c828005fd24a55b519e68b02  -

