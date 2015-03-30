Set up
  $ . $TESTDIR/setup.sh

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
  $ sed -n 6,110864p $STDDIR/concordant_subset_2015_03_28.sam > $OUTDIR/tmp2
  $ sort $OUTDIR/tmp2 > $OUTDIR/tmp22
  $ diff $OUTDIR/tmp11 $OUTDIR/tmp22
  $ rm -rf $OUTDIR/tmp1 $OUTDIR/tmp2 $OUTDIR/tmp11 $OUTDIR/tmp22
#concordant_subset_2014_05_28.sam  --> changelist 135254, use MAX_BAND_SIZE to contrain GuidedAlign
#concordant_subset_2014_08_21.sam  --> changelist 138516, added YS, YE, ZM tags. 
#concordant_subset_2014_08_28.sam  --> changelist 139176, update SAM MD5 
#concordant_subset_2014_09_12.sam  --> changelist 140410, changed the default value of '-concordantTemplate' from 'longestsubread' to 'typicalsubread'
#concordant_subset_2014_09_17.sam  --> changelist 140573, changed SDPFragment LessThan to make sure blasr compiled with gcc 4.4 and 4.8 can produce identical results. 
#concordant_subset_2014_10_16.sam  --> changelist 141378, changed the default value of '-concordantTemplate' from 'typicalsubread' to 'mediansubread'
#concordant_subset_2015_03_01.sam  --> changelist 146599, reads from the same movie should have unique readGroupId
#concordant_subset_2015_03_28.sam  --> changelist 148101, 148080 updated read group id, 148100 updated TLEN

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
