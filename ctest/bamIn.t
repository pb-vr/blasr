Set up
  $ . $TESTDIR/setup.sh

Test using bam as input 
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/lambda_ref.fasta -m 4 --out $OUTDIR/tiny_bam_in.m4 
  [INFO]* (glob)
  [INFO]* (glob)

Check whether blasr produces identical results taking fasta sequences of the bam as input 
  $ $EXEC $DATDIR/test_bam/tiny_fasta.fofn $DATDIR/lambda_ref.fasta -m 4 --out $OUTDIR/tiny_fasta_in.m4
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/tiny_bam_in.m4 $OUTDIR/tiny_fasta_in.m4

Test bam in, bam out
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/lambda_ref.fasta --bam --out $OUTDIR/tiny_bam_in.bam --clipping subread
  [INFO]* (glob)
  [INFO]* (glob)

TODO: test --concordant, when pbbam API to query over ZMWs is available.
TODO: test bam with ccs reads
