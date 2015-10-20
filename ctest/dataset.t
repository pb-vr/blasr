Set up
  $ . $TESTDIR/setup.sh

Must copy all data to the current dir before bug 28912 is solved.
  $ cp $DATDIR/test_dataset/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam . 
  $ cp $DATDIR/test_dataset/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam.pbi . 
  $ cp $DATDIR/test_dataset/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.2.subreads.bam . 
  $ cp $DATDIR/test_dataset/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.2.subreads.bam.pbi . 
  $ cp $DATDIR/test_dataset/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.3.subreads.bam . 
  $ cp $DATDIR/test_dataset/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.3.subreads.bam.pbi . 

Test dataset.xml as input
  $ $EXEC $DATDIR/test_dataset/chunking.subreadset.xml $DATDIR/ecoli_reference.fasta -m 4 -out $OUTDIR/chunking.m4 -bestn 1 && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0
Test filters in dataset.xml is respected.
  $ cat $OUTDIR/chunking.m4 | wc -l
  9

Test dataset.xml -bam output
  $ $EXEC $DATDIR/test_dataset/chunking.subreadset.xml $DATDIR/ecoli_reference.fasta -bam -out $OUTDIR/chunking.bam  && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0

Test dataset.xml -concordant
  $ $EXEC $DATDIR/test_dataset/chunking.subreadset.xml $DATDIR/ecoli_reference.fasta -bam -out $OUTDIR/chunking.concordant.bam -concordant && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0

Test dataset with no filters (to make sure that an empty filter does not discard all bam records.)
  $ $EXEC $DATDIR/test_dataset/nofilter.subreadset.xml $DATDIR/ecoli_reference.fasta -bam -out $OUTDIR/nofilter.bam -concordant -bestn 1 && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0

  $ $SAMTOOLS view $OUTDIR/nofilter.bam|wc -l
  135
