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
  [INFO]* (glob)
  [INFO]* (glob)
  $ tail -n+5 $OUTDIR/ecoli.sam | sort | cut -f 1-11| md5sum
  d2e2b6cfe710b7b6a065e73c3244b0b9  -

#  $ tail -n+5 $STDDIR/ecoli.sam | sort | cut -f 1-11| md5sum
#  d2e2b6cfe710b7b6a065e73c3244b0b9  -

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
  $ diff $OUTDIR/read.m1 $STDDIR/read.m1

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

Test blasr with *.fofn input
  $ rm -rf $OUTDIR/lambda_bax.m4 
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out lambda_bax_tmp.m4 -nproc 15 -minMatch 14
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort lambda_bax_tmp.m4 > $OUTDIR/lambda_bax.m4
  $ diff $OUTDIR/lambda_bax.m4 $STDDIR/lambda_bax.m4

Test blasr with -noSplitSubreads 
  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads.m4 $OUTDIR/lambda_bax_noSplitSubreads.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -noSplitSubreads -m 4 -out lambda_bax_noSplitSubreads_tmp.m4 -nproc 15
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort lambda_bax_noSplitSubreads_tmp.m4 > $OUTDIR/lambda_bax_noSplitSubreads.m4
  $ diff $OUTDIR/lambda_bax_noSplitSubreads.m4 $STDDIR/lambda_bax_noSplitSubreads.m4

Test alignment score
  $ rm -rf $OUTDIR/testscore.m0
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -minReadLength 1 -m 0 -out $OUTDIR/testscore.m0
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/testscore.m0 $STDDIR/testscore.m0

Test affineAlign
  $ rm -rf $OUTDIR/affineAlign.m0
  $ $EXEC $DATDIR/affineAlign.fofn $DATDIR/substr_with_ins.fasta -m 0 -out $OUTDIR/affineAlign.m0  -affineAlign  -readIndex 493 -insertion 100 -deletion 100
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/affineAlign.m0 $STDDIR/affineAlign.m0

  $ rm -rf $OUTDIR/ecoli_affine.m0
  $ $EXEC $DATDIR/ecoli_affine.fasta $DATDIR/ecoli_reference.fasta -m 0 -out $OUTDIR/ecoli_affine.m0 -affineAlign -insertion 100 -deletion 100
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/ecoli_affine.m0 $STDDIR/ecoli_affine.m0


Test -holeNumbers
  $ rm -f $OUTDIR/holeNumbers.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out $OUTDIR/holeNumbers.m4 -holeNumbers 14798,55000-55100 -nproc 8
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/holeNumbers.m4 | md5sum
  21fd37b14b85ef7dda332ea10edc524a  -

Test Sam out nm tag
  $ rm -rf $OUTDIR/read.sam
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -sam -out $OUTDIR/read.sam
  [INFO]* (glob)
  [INFO]* (glob)
  $ tail -n+5 $OUTDIR/read.sam |cut -f 18 
  NM:i:2
  NM:i:3
  NM:i:2
  NM:i:4

Test -useccsall with bestn = 1
  $ $EXEC $DATDIR/ccstest.fofn $DATDIR/ccstest_ref.fasta -bestn 1 -useccsall -sam -out $OUTDIR/useccsall.sam -holeNumbers 76772
  [INFO]* (glob)
  [INFO]* (glob)
  $ tail -n+9 $OUTDIR/useccsall.sam | md5sum 
  8f9cba19d956a140b359e89db714bbf4  -


Test -concordant
  $ rm -rf $OUTDIR/concordant.sam
  $ $EXEC $DATDIR/ecoli_lp.fofn $DATDIR/ecoli_reference.fasta -concordant -sam -out $OUTDIR/concordant.sam -nproc 8
  [INFO]* (glob)
  [INFO]* (glob)
  $ sed -n 6,110864p $OUTDIR/concordant.sam > $OUTDIR/tmp1 
  $ sort $OUTDIR/tmp1 > $OUTDIR/tmp11
  $ sed -n 6,110864p $STDDIR/concordant.sam > $OUTDIR/tmp2
  $ sort $OUTDIR/tmp2 > $OUTDIR/tmp22
  $ diff $OUTDIR/tmp11 $OUTDIR/tmp22

  $ rm -rf $OUTDIR/tmp1 $OUTDIR/tmp2 $OUTDIR/tmp11 $OUTDIR/tmp22

Test -concordant 
  $ rm -f $OUTDIR/concordant2.samtom4 $OUTDIR/concordant2.sam $OUTDIR/not_concordant2.m4
  $ FOFN=$DATDIR/concordant.fofn
  $ REF=$DATDIR/lambda_ref.fasta
  $ $EXEC $FOFN $REF -concordant -sam -out $OUTDIR/concordant2.sam -holeNumbers 4405
  [INFO]* (glob)
  [INFO]* (glob)
  $ $EXEC $FOFN $REF -m 4 -out $OUTDIR/not_concordant2.m4 -holeNumbers 4405
  [INFO]* (glob)
  [INFO]* (glob)
  $ $TESTDIR/../pbihdfutils/bin/samtom4 $OUTDIR/concordant2.sam $REF $OUTDIR/concordant2.samtom4 
  $ diff $OUTDIR/not_concordant2.m4 $OUTDIR/concordant2.samtom4


Test using *.ccs.h5 as input
# The results should be exactly the same as 
# blasr ccsasinput_bas.fofn $DATDIR/ccsasinput.fasta -m 4 -out tmp.m4 -useccsdenovo
  $ rm -rf $OUTDIR/ccsasinput.m4
  $ $EXEC $DATDIR/ccsasinput.fofn $DATDIR/ccsasinput.fasta -m 4 -out $OUTDIR/ccsasinput.m4
  [INFO]* (glob)
  [INFO]* (glob)
  $ cat $OUTDIR/ccsasinput.m4 | md5sum
  55295b4304c1cd1e79edb810bb048a4c  -

Test -useccsall with Large genome.
  $ BASFILE=/mnt/data3/vol53/2450530/0014/Analysis_Results/m130507_052228_42161_c100519212550000001823079909281305_s1_p0.3.bax.h5
  $ REFDIR=/mnt/secondary/Smrtpipe/repository/hg19_M_sorted/sequence
  $ REFFA=$REFDIR/hg19_M_sorted.fasta
  $ REFSA=$REFDIR/hg19_M_sorted.fasta.sa
  $ OUTFILE=$OUTDIR/intflow.m4
  $ $EXEC $BASFILE $REFFA -out $OUTFILE -m 4 -sa $REFSA -holeNumbers 109020
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $STDDIR/intflow.m4 $OUTFILE


