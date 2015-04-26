Set up
  $ . $TESTDIR/setup.sh

Test generating bam output 
Input is bax.h5
  $ $EXEC $DATDIR/test_bam/tiny_bax.fofn $DATDIR/lambda_ref.fasta -bam -out $OUTDIR/tiny_bax_in.bam -clipping soft
  [INFO]* (glob)
  [INFO]* (glob)

  $ $SAMTOOLS view -h $OUTDIR/tiny_bax_in.bam -o $OUTDIR/tiny_bax_in.bam.sam

  $ $EXEC $DATDIR/test_bam/tiny_bax.fofn $DATDIR/lambda_ref.fasta -sam -out $OUTDIR/tiny_bax_in.sam -clipping soft
  [INFO]* (glob)
  [INFO]* (glob)

Check whether sam and bam produces identical alignments
  $ sed -n '6,$p' $OUTDIR/tiny_bax_in.bam.sam | cut -f 2-11 > $TMP1.bax_in
  $ sed -n '6,$p' $OUTDIR/tiny_bax_in.sam | cut -f 2-11 > $TMP2.bax_in
  $ diff $TMP1.bax_in $TMP2.bax_in

Compare with stdout
  $ sed -n '6,$p' $STDDIR/$UPDATEDATE/tiny_bax_in.bam.sam | cut -f 2-11 > $TMP2.bax_in
  $ diff $TMP1.bax_in $TMP2.bax_in

Input is fasta, compare with stdout
  $ $EXEC $DATDIR/test_bam/tiny_fasta.fofn $DATDIR/lambda_ref.fasta -bam -out $OUTDIR/tiny_fasta_in.bam -clipping soft
  [INFO]* (glob)
  [INFO]* (glob)

  $ $SAMTOOLS view -h $OUTDIR/tiny_fasta_in.bam -o $OUTDIR/tiny_fasta_in.bam.sam
  $ sed -n '6,$p' $OUTDIR/tiny_fasta_in.bam.sam > $TMP1.fasta_in
  $ sed -n '6,$p' $STDDIR/$UPDATEDATE/tiny_fasta_in.bam.sam > $TMP2.fasta_in
  $ diff $TMP1.fasta_in $TMP2.fasta_in

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

