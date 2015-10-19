Set up
  $ . $TESTDIR/setup.sh

Test using bam as input, use -concordant
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/bamConcordantRef.fasta -bam -concordant -bestn 1 -out $OUTDIR/bamConcordant.bam
  [INFO]* (glob)
  [INFO]* (glob)

Check whether sam out and bam out have identical alignments, not checking qvs
  $ $SAMTOOLS view $OUTDIR/bamConcordant.bam |cut -f 4 
  8067
  8051
  730
  690
  690
  690
  691
  690
  690
  693
  690
  690
  690
  697
  695
  690
