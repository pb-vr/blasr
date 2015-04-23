Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: pls2fasta.
  $ EXEC=$TESTDIR/../pls2fasta

Test pls2fasta
Condition: the order of region tables do not match the order bax.h5 files. 
  $ $EXEC $DATDIR/test_pls2fasta.fofn $OUTDIR/test_pls2fasta.fa -regionTable $DATDIR/test_pls2fasta_rgn.fofn -trimByRegion
  [INFO] * [pls2fasta] started. (glob)
  [INFO] * [pls2fasta] ended. (glob)


Test pls2fasta output fasta
  $ $EXEC $DATDIR/ecoli_lp.fofn $OUTDIR/test_pls2fasta_ecoli.fa -trimByRegion
  [INFO] * [pls2fasta] started. (glob)
  [INFO] * [pls2fasta] ended. (glob)
  $ echo $?
  0
  $ diff $OUTDIR/test_pls2fasta_ecoli.fa $STDDIR/test_pls2fasta_ecoli.fa

Test pls2fasta output fastq
  $ $EXEC $DATDIR/ecoli_lp.fofn $OUTDIR/test_pls2fasta_ecoli.fq -trimByRegion -fastq
  [INFO] * [pls2fasta] started. (glob)
  [INFO] * [pls2fasta] ended. (glob)
  $ echo $?
  0
  $ diff $OUTDIR/test_pls2fasta_ecoli.fq $STDDIR/test_pls2fasta_ecoli.fq
