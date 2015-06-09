Set up
  $ . $TESTDIR/setup.sh

Test fasta input, bam out
  $ $EXEC $DATDIR/tiny_fasta.fofn $DATDIR/lambda_ref.fasta -bam -out $OUTDIR/fasta_in_bam_out.bam
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0

  $ $EXEC $DATDIR/tiny_fasta.fofn $DATDIR/lambda_ref.fasta -bam -out $OUTDIR/fasta_in_bam_out2.bam -noSplitSubreads
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
