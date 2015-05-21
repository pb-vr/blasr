Set up
  $ . $TESTDIR/setup.sh

  $ NAME=test_filtercriteria
  $ DATDIR=$DATDIR/$NAME
  $ OUTDIR=$OUTDIR/$NAME
  $ STDDIR=$STDDIR/$NAME
  $ mkdir -p $OUTDIR

Test -minPctSimilarity
  $ I=$DATDIR/tiny_bam.fofn
  $ R=$DATDIR/lambdaNEB.fa
  $ O=$OUTDIR/min_pct_similarity_90.m4

  $ $EXEC $I $R -out $O -m 4 -minPctSimilarity 90
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
  $ awk '$4 < 90 {print}' |wc -l |cut -f 1 -d ' ' 
  0

  $ O=$OUTDIR/min_aln_len_1000.m4
  $ $EXEC $I $R -out $O -m 4 -minAlnLength 1000
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
  $ wc -l $O |cut -f 1 -d ' ' 
  12
