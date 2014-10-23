Set up directories
  $ CURDIR=$TESTDIR
  $ REMOTEDIR=/mnt/secondary-siv/testdata/BlasrTestData/ctest
  $ DATDIR=$REMOTEDIR/data
  $ OUTDIR=$CURDIR/out
  $ STDDIR=$REMOTEDIR/stdout

Set up the executable: blasr.
  $ EXEC=$TESTDIR/../blasr

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

#string perforceVersionString("$Change: 139404 $");

Test blasr on ecoli.
Test blasr with -sam

# The following job takes a very long time to finish, let us use a subset of reads instead
#See $STDOUT/ecoli_v1.4.sam for 1.4 output.
# $STDOUT/ecoli_2014_03_28.sam for bug before mapQV for affineAlign/align without QV is fixed.

#  $ rm -rf $OUTDIR/ecoli.sam
#  $ $EXEC $DATDIR/ecoli.fasta $DATDIR/ecoli_reference.fasta -sam -out $OUTDIR/ecoli.sam -nproc 15
#  [INFO]* (glob)
#  [INFO]* (glob)
#
#  $ sed -n '5,$ p' $OUTDIR/ecoli.sam | sort | cut -f 1-11 > $TMP1
#  $ sed -n '5,$ p' $STDDIR/ecoli_2014_05_04.sam | sort | cut -f 1-11 > $TMP2
#  $ diff $TMP1 $TMP2
#  $ rm $TMP1 $TMP2

  $ rm -rf $OUTDIR/ecoli_subset.sam
  $ $EXEC $DATDIR/ecoli_subset.fasta $DATDIR/ecoli_reference.fasta -sam -out $OUTDIR/ecoli_subset.sam -nproc 15
  [INFO]* (glob)
  [INFO]* (glob)

  $ sed -n '5,$ p' $OUTDIR/ecoli_subset.sam | sort | cut -f 1-11 > $TMP1
  $ sed -n '5,$ p' $STDDIR/ecoli_subset_2014_06_10.sam | sort | cut -f 1-11 > $TMP2
  $ diff $TMP1 $TMP2
  $ rm $TMP1 $TMP2


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


Test blasr with *.fofn input
#  $ rm -rf $OUTDIR/lambda_bax.m4 
#  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out lambda_bax_tmp.m4 -nproc 15 -minMatch 14
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ sort lambda_bax_tmp.m4 > $OUTDIR/lambda_bax.m4
#  $ diff $OUTDIR/lambda_bax.m4 $STDDIR/lambda_bax.m4
# This test takes a long time, use a subset instad. 

  $ rm -rf $OUTDIR/lambda_bax_subset.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out $OUTDIR/lambda_bax_tmp_subset.m4 -nproc 15 -minMatch 14 -holeNumbers 1-1000 -sa $DATDIR/lambda_ref.sa
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/lambda_bax_tmp_subset.m4 > $OUTDIR/lambda_bax_subset.m4
  $ diff $OUTDIR/lambda_bax_subset.m4 $STDDIR/lambda_bax_subset.m4


Test blasr with -noSplitSubreads 
#  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads.m4 
#  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -noSplitSubreads -m 4 -out lambda_bax_noSplitSubreads_tmp.m4 -nproc 15
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ sort lambda_bax_noSplitSubreads_tmp.m4 > $OUTDIR/lambda_bax_noSplitSubreads.m4
#  $ diff $OUTDIR/lambda_bax_noSplitSubreads.m4 $STDDIR/lambda_bax_noSplitSubreads.m4
# This test takes a long time, use a subset instad. 

  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads_subset.m4 
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -noSplitSubreads -m 4 -out $OUTDIR/lambda_bax_noSplitSubreads_tmp_subset.m4 -nproc 15 -holeNumbers 1-1000 -sa $DATDIR/lambda_ref.sa
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/lambda_bax_noSplitSubreads_tmp_subset.m4 > $OUTDIR/lambda_bax_noSplitSubreads_subset.m4
  $ diff $OUTDIR/lambda_bax_noSplitSubreads_subset.m4 $STDDIR/lambda_bax_noSplitSubreads_subset.m4


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
  $ diff $OUTDIR/affineAlign.m0 $STDDIR/affineAlign_2014_06_10.m0


  $ rm -rf $OUTDIR/ecoli_affine.m0
  $ $EXEC $DATDIR/ecoli_affine.fasta $DATDIR/ecoli_reference.fasta -m 0 -out $OUTDIR/ecoli_affine.m0 -affineAlign -insertion 100 -deletion 100
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/ecoli_affine.m0 $STDDIR/ecoli_affine_2014_06_10.m0
# Note that MapQV for -affineAlign has been fixed in 2014 04 18, bug 24363 

Test -holeNumbers
  $ rm -f $OUTDIR/holeNumbers.m4
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta -m 4 -out $OUTDIR/holeNumbers.m4 -holeNumbers 14798,55000-55100 -nproc 8
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/holeNumbers.m4 > $TMP1
  $ sort $STDDIR/holeNumbers_2014_05_29.m4 > $TMP2
  $ diff $TMP1 $TMP2
  $ rm $TMP1 $TMP2

Test Sam out nm tag
  $ rm -rf $OUTDIR/read.sam
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -sam -out $OUTDIR/read.sam
  [INFO]* (glob)
  [INFO]* (glob)
  $ tail -n+5 $OUTDIR/read.sam |cut -f 21 
  NM:i:2
  NM:i:3
  NM:i:2
  NM:i:4

Test -useccsall with bestn = 1
  $ $EXEC $DATDIR/ccstest.fofn $DATDIR/ccstest_ref.fasta -bestn 1 -useccsall -sam -out $OUTDIR/useccsall.sam -holeNumbers 76772
  [INFO]* (glob)
  [INFO]* (glob)
  $ sed -n '9,$ p' $OUTDIR/useccsall.sam > $TMP1
  $ sed -n '9,$ p' $STDDIR/useccsall_2014_08_28.sam > $TMP2
  $ diff $TMP1 $TMP2
  $ rm $TMP1 $TMP2
# useccsall_2014_08_21.sam --> change 138516: added YS, YE, ZM tags
# useccsall_2014_08_28.sam --> change 139176: Update SAM MD5 

Test -concordant
#  $ rm -rf $OUTDIR/concordant.sam
#  $ $EXEC $DATDIR/ecoli_lp.fofn $DATDIR/ecoli_reference.fasta -concordant -sam -out $OUTDIR/concordant.sam -nproc 8
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ sed -n 6,110864p $OUTDIR/concordant.sam > $OUTDIR/tmp1 
#  $ sort $OUTDIR/tmp1 > $OUTDIR/tmp11
#  $ sed -n 6,110864p $STDDIR/concordant.sam > $OUTDIR/tmp2
#  $ sort $OUTDIR/tmp2 > $OUTDIR/tmp22
#  $ diff $OUTDIR/tmp11 $OUTDIR/tmp22
#  $ rm -rf $OUTDIR/tmp1 $OUTDIR/tmp2 $OUTDIR/tmp11 $OUTDIR/tmp22
# This test takes a long time, use a subset instad. 

  $ rm -rf $OUTDIR/concordant_subset.sam
  $ $EXEC $DATDIR/ecoli_lp.fofn $DATDIR/ecoli_reference.fasta -concordant -sam -out $OUTDIR/concordant_subset.sam -nproc 12 -holeNumbers 1-10000 -sa $DATDIR/ecoli_reference.sa
  [INFO]* (glob)
  [INFO]* (glob)
  $ sed -n 6,110864p $OUTDIR/concordant_subset.sam > $OUTDIR/tmp1 
  $ sort $OUTDIR/tmp1 > $OUTDIR/tmp11
  $ sed -n 6,110864p $STDDIR/concordant_subset_2014_10_16.sam > $OUTDIR/tmp2
  $ sort $OUTDIR/tmp2 > $OUTDIR/tmp22
  $ diff $OUTDIR/tmp11 $OUTDIR/tmp22
  $ rm -rf $OUTDIR/tmp1 $OUTDIR/tmp2 $OUTDIR/tmp11 $OUTDIR/tmp22
#concordant_subset_2014_05_28.sam  --> changelist 135254, use MAX_BAND_SIZE to contrain GuidedAlign
#concordant_subset_2014_08_21.sam  --> changelist 138516, added YS, YE, ZM tags. 
#concordant_subset_2014_08_28.sam  --> changelist 139176, update SAM MD5 
#concordant_subset_2014_09_12.sam  --> changelist 140410, changed the default value of '-concordantTemplate' from 'longestsubread' to 'typicalsubread'
#concordant_subset_2014_09_17.sam  --> changelist 140573, changed SDPFragment LessThan to make sure blasr compiled with gcc 4.4 and 4.8 can produce identical results. 
#concordant_subset_2014_10_16.sam  --> changelist 141378, changed the default value of '-concordantTemplate' from 'typicalsubread' to 'mediansubread'

#Test -concordant, case 2
#Avoid to use samtom4 in blasr cram tests, since blasr has been seprated 
# from other executables, and it is good to make blasr cram test not depend
# on other tools, but only focus on blasr specific issues.
#  $ rm -f $OUTDIR/concordant2.samtom4 $OUTDIR/concordant2.sam $OUTDIR/not_concordant2.m4
#  $ FOFN=$DATDIR/concordant.fofn
#  $ REF=$DATDIR/lambda_ref.fasta
#  $ $EXEC $FOFN $REF -concordant -sam -out $OUTDIR/concordant2.sam -holeNumbers 4405
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ $EXEC $FOFN $REF -m 4 -out $OUTDIR/not_concordant2.m4 -holeNumbers 4405
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ samtom4 $OUTDIR/concordant2.sam $REF $OUTDIR/concordant2.samtom4 
#  $ diff $OUTDIR/not_concordant2.m4 $OUTDIR/concordant2.samtom4

Test -concordant, case 2
  $ rm -f $OUTDIR/concordant2.m4 $OUTDIR/not_concordant2.m4
  $ FOFN=$DATDIR/concordant.fofn
  $ REF=$DATDIR/lambda_ref.fasta
  $ $EXEC $FOFN $REF -concordant -out $OUTDIR/concordant2.m4 -holeNumbers 4405
  [INFO]* (glob)
  [INFO]* (glob)
  $ $EXEC $FOFN $REF -m 4 -out $OUTDIR/not_concordant2.m4 -holeNumbers 4405
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/concordant2.m4 $STDDIR/concordant2.m4 
  $ diff $OUTDIR/not_concordant2.m4 $STDDIR/not_concordant2.m4

Test -concordant FMR1 case (the 'typical subread' is selected as template for concordant mapping)
  $ FOFN=$DATDIR/FMR1_concordant.fofn
  $ REF=$DATDIR/FMR1_130CGG.fasta
  $ $EXEC $FOFN $REF -concordant -out $OUTDIR/FMR1_zmw_37927.m4 -m 4 -holeNumbers 37927
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/FMR1_zmw_37927.m4 $STDDIR/FMR1_zmw_37927.m4

Test using *.ccs.h5 as input
# The results should be exactly the same as 
# blasr $DATDIR/ccsasinput_bas.fofn $DATDIR/ccsasinput.fasta -m 4 -out tmp.m4 -useccsdenovo
  $ rm -rf $OUTDIR/ccsasinput.m4
  $ $EXEC $DATDIR/ccsasinput.fofn $DATDIR/ccsasinput.fasta -m 4 -out $OUTDIR/ccsasinput.m4
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/ccsasinput.m4 $STDDIR/ccsasinput_2014_06_10.m4

Test -useccsall with Large genome.
  $ BASFILE=/mnt/data3/vol53/2450530/0014/Analysis_Results/m130507_052228_42161_c100519212550000001823079909281305_s1_p0.3.bax.h5
  $ REFDIR=/mnt/secondary/Smrtpipe/repository/hg19_M_sorted/sequence
  $ REFFA=$REFDIR/hg19_M_sorted.fasta
  $ REFSA=$REFDIR/hg19_M_sorted.fasta.sa
  $ OUTFILE=$OUTDIR/intflow.m4
  $ $EXEC $BASFILE $REFFA -out $OUTFILE -m 4 -sa $REFSA -holeNumbers 109020
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $STDDIR/intflow_2014_06_10.m4 $OUTFILE


Test input.fofn containing a new bas.h5 file. Note that the new bas.h5 file does not 
contain any /PulseData, instead contains /MultiPart/Parts.
  $ rm -f $TMP1
  $ BASFILE=/mnt/data3/vol53/2450598/0001/Analysis_Results/m130812_185809_42141_c100533960310000001823079711101380_s1_p0.bas.h5
  $ REFFA=/mnt/secondary/Smrtpipe/repository/Ecoli_BL21_O26/sequence/Ecoli_BL21_O26.fasta
  $ $EXEC $BASFILE $REFFA -holeNumbers 1-100 -out $TMP1
  [INFO] * [blasr] started. (glob)
  [INFO] * [blasr] ended. (glob)
  $ echo $?
  0

Test -fastMaxInterval.
  $ rm -f $TMP1
  $ BASFILE=/mnt/data3/vol53/2450598/0001/Analysis_Results/m130812_185809_42141_c100533960310000001823079711101380_s1_p0.bas.h5
  $ REFFA=/mnt/secondary/Smrtpipe/repository/Ecoli_BL21_O26/sequence/Ecoli_BL21_O26.fasta
  $ $EXEC $BASFILE $REFFA -holeNumbers 1-100 -out $TMP1 -fastMaxInterval
  [INFO] * [blasr] started. (glob)
  [INFO] * [blasr] ended. (glob)
  $ echo $?
  0

Test -aggressiveIntervalCut.
  $ rm -f $TMP1
  $ BASFILE=/mnt/data3/vol53/2450598/0001/Analysis_Results/m130812_185809_42141_c100533960310000001823079711101380_s1_p0.bas.h5
  $ REFFA=/mnt/secondary/Smrtpipe/repository/Ecoli_BL21_O26/sequence/Ecoli_BL21_O26.fasta
  $ $EXEC $BASFILE $REFFA -holeNumbers 1-100 -out $TMP1 -aggressiveIntervalCut
  [INFO] * [blasr] started. (glob)
  [INFO] * [blasr] ended. (glob)
  $ echo $?
  0

bug_25328, unrolled resequencing test 
  $ INFA=$DATDIR/bug_25328_zmw_38131.fasta
  $ REF=$DATDIR/All4mers_circular_72x_l50256.fasta
  $ OUTFA=$OUTDIR/bug_25328.m4
  $ $EXEC $INFA $REF -bestn 1 -nCandidates 1 -forwardOnly -maxMatch 14 -m 4 -out $OUTFA
  [INFO]* (glob)
  [INFO]* (glob)

  $ awk '$7-$6 >= 15000' $OUTFA |wc -l
  1

bug_25766, added an option -minRawSubreadScore
  $ BASFILE=$DATDIR/lambda_bax.fofn
  $ REF=$DATDIR/lambda_ref.fasta
  $ $EXEC $BASFILE $REF -out $TMP1 -minRawSubreadScore 700 -nproc 18
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
  $ sort $TMP1 > $TMP2
  $ diff $TMP2 $STDDIR/bug_25766.m4

bug 25741 If input bas.h5 does not contain mergeQV, blasr with -printSAMQV, -nproc>1 should not write garbage 'mq' values to output.
  $ $EXEC $DATDIR/bas_wo_mergeQV.fofn $DATDIR/lambda_ref.fasta -printSAMQV -sam -clipping subread -out $OUTDIR/out_printSAMQV.sam -nproc 12 
  [INFO]* (glob)
  [INFO]* (glob)
  $ grep 'mq' $OUTDIR/out_printSAMQV.sam |wc -l
  1

