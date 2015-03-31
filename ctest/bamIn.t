Set up
  $ . $TESTDIR/setup.sh

Test using bam as input 
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/lambda_ref.fasta -m 4 -out $OUTDIR/tiny_bam_in.m4 
  [INFO]* (glob)
  [INFO]* (glob)

Check whether blasr produces identical results taking fasta sequences of the bam as input 
  $ $EXEC $DATDIR/test_bam/tiny_fasta.fofn $DATDIR/lambda_ref.fasta -m 4 -out $OUTDIR/tiny_fasta_in.m4
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/tiny_bam_in.m4 $OUTDIR/tiny_fasta_in.m4

Test bam in, sam out
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/lambda_ref.fasta -sam -out $OUTDIR/tiny_bam_in.sam -printSAMQV -clipping subread
  [INFO]* (glob)
  [INFO]* (glob)

Test bam in, bam out
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/lambda_ref.fasta -bam -out $OUTDIR/tiny_bam_in.bam -clipping subread
  [INFO]* (glob)
  [INFO]* (glob)

Check whether sam out and bam out have identical alignments, not checking qvs
  $ $SAMTOOLS view -h $OUTDIR/tiny_bam_in.bam -o $OUTDIR/tiny_bam_in.bam.sam
  $ cut -f 2-11 $OUTDIR/tiny_bam_in.bam.sam |sed -n '6,$p' > $TMP1.aln
  $ cut -f 2-11 $OUTDIR/tiny_bam_in.sam |sed -n '6,$p' > $TMP2.aln
  $ diff $TMP1.aln $TMP2.aln

Check whether sam out and bam out have identical read groups @RG 
  $ awk '/^@RG/' $OUTDIR/tiny_bam_in.bam.sam > $TMP1.rg
  $ awk '/^@RG/' $OUTDIR/tiny_bam_in.sam > $TMP2.rg
  $ diff $TMP1.rg $TMP2.rg

Compare iq produced with stdout
  $ sed -n '6,$p' $OUTDIR/tiny_bam_in.bam.sam | awk '{gsub(/\t/,"\n");}1' | awk '/^iq:Z:/' > $TMP1.iq
  $ sed -n '6,$p' $STDDIR/tiny_bam_in.bam.sam | awk '{gsub(/\t/,"\n");}1' | awk '/^iq:Z:/' > $TMP2.iq
  $ diff $TMP1.iq $TMP2.iq

TODO:Check whether sam out and bam out have identical insertion qvs
Currently QVs in bam are in 'native' orientation, and QVs in sam are in 'genomic' orientation. This needs to be fixed.
$ sed -n '6,$p' $OUTDIR/tiny_bam_in.sam | awk '{gsub(/\t/,"\n");}1' | awk '/^iq:Z:/' > $TMP2.iq

Test with multiple nproc
  $ $EXEC $DATDIR/test_bam/two_bam.fofn $DATDIR/lambda_ref.fasta -bam -nproc 15 -out $OUTDIR/two_bam_in.bam 
  [INFO]* (glob)
  [INFO]* (glob)
  $ $SAMTOOLS view -h $OUTDIR/two_bam_in.bam -o $OUTDIR/two_bam_in.bam.sam

TODO: test -concordant, when pbbam API to query over ZMWs is available.
TODO: test bam with ccs reads
