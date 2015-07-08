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

Test if bam cigar strings are correct
  $ head -2 $TMP1.bam_in_soft |cut -f 6
  25=1I28=1I41=1I5=1D6=1X12=1I15=1I2=1I16=1D10=1I11=1I74=1D12=1D7=3I4=1I6=1D1=2D14=1D16=1I8=1D4=1D5=1D20=1I3=1I10=1I37=1I13=1I25=1I15=1I7=1I11=1I3=2I1=1I16=1I6=1I8=1I11=1X1=1I5=1I56=1I17=
  28=1D7=1I1=1I9=2I12=1I3=1D13=1I15=1I2=1X49=1I19=1I14=1I5=1D17=1D20=1D86=1I21=1I9=1I24=1I6=1I1=1I2=1D11=1D4=1D3=1D31=1D6=1I6=1I9=1I57=2I24=1I26=1I8=1I43=1S
