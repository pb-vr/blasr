Set up
  $ . $TESTDIR/setup.sh

Test generating bam output 

Input is bam, clipping=soft and subread should produce identical results
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/lambda_ref.fasta -bam -out $OUTDIR/tiny_bam_in_soft.bam -clipping soft
  [INFO]* (glob)
  [INFO]* (glob)

  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/lambda_ref.fasta -bam -out $OUTDIR/tiny_bam_in_subread.bam -clipping subread 
  [INFO]* (glob)
  [INFO]* (glob)

  $ $SAMTOOLS view $OUTDIR/tiny_bam_in_soft.bam | sed -n '6,$p' > $TMP1.bam_in_soft
  $ $SAMTOOLS view $OUTDIR/tiny_bam_in_subread.bam | sed -n '6,$p' > $TMP2.bam_in_subread
  $ diff $TMP1.bam_in_soft $TMP2.bam_in_subread

