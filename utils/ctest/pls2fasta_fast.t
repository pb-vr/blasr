Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: pls2fasta.
  $ EXEC=$TESTDIR/../pls2fasta

Test pls2fasta output fasta
  $ $EXEC $DATDIR/ecoli_lp.fofn $OUTDIR/test_pls2fasta_ecoli.fa -trimByRegion
  [INFO] * [pls2fasta] started. (glob)
  [INFO] * [pls2fasta] ended. (glob)
  $ echo $?
  0
  $ diff $OUTDIR/test_pls2fasta_ecoli.fa $STDDIR/test_pls2fasta_ecoli.fa


