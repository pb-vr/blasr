Set up
  $ . $TESTDIR/setup.sh

  $ NAME=test_hitpolicy
  $ DATDIR=$DATDIR/$NAME
  $ OUTDIR=$OUTDIR/$NAME
  $ STDDIR=$STDDIR/$NAME
  $ mkdir -p $OUTDIR

  $ I=$DATDIR/tiny_bam.fofn
  $ R=$DATDIR/test_hitpolicy_target.fa
  $ O=$OUTDIR/hitpolicy_all.m4
  $ X=$STDDIR/hitpolicy_all.m4

Test hitpolicy all 
  $ $EXEC $I $R -out $O -m 4 -hitPolicy all
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
  $ wc -l $O | cut -f 1 -d ' '
  683

Test hitpolicy allbest 
  $ O=$OUTDIR/hitpolicy_allbest.m4
  $ X=$STDDIR/hitpolicy_allbest.m4
  $ $EXEC $I $R -out $O -m 4 -hitPolicy allbest && sort $O > $TMP1 && mv $TMP1 $O
  [INFO]* (glob)
  [INFO]* (glob)
  $ echo $?
  0
  $ sort $O > $TMP1 && mv $TMP1 $O
  $ diff $O $X && echo $?
  0

Test hitpolicy random
  $ O=$OUTDIR/hitpolicy_random.m4
  $ X=$STDDIR/hitpolicy_random.m4
  $ $EXEC $I $R -out $O -m 4 -hitPolicy random -randomSeed 1
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $O > $TMP1 && mv $TMP1 $O
  $ diff $O $X && echo $?
  0

Test hitpolicy randombest bam inputs, nproc > 1, fixed seed
  $ O=$OUTDIR/hitpolicy_randombest_bam_in.m4
  $ X=$STDDIR/hitpolicy_randombest_bam_in.m4
  $ $EXEC $I $R -out $O -m 4 -hitPolicy randombest -randomSeed 1 -nproc 10
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $O > $TMP1 && mv $TMP1 $O
  $ diff $O $X && echo $?
  0

Test hitpolicy randombest bax inputs, nproc > 1, fixed seed
  $ I=$DATDIR/tiny_bax.fofn
  $ O=$OUTDIR/hitpolicy_randombest_bax_in.m4
  $ X=$STDDIR/hitpolicy_randombest_bax_in.m4
  $ $EXEC $I $R -out $O -m 4 -hitPolicy randombest -randomSeed 1 -nproc 10
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $O > $TMP1 && mv $TMP1 $O
  $ diff $O $X && echo $?
  0

Test hitpolicy randombest fasta inputs, nproc > 1, fixed seed
  $ I=$DATDIR/tiny_fasta.fofn
  $ O=$OUTDIR/hitpolicy_randombest_fasta_in.m4
  $ X=$STDDIR/hitpolicy_randombest_fasta_in.m4
  $ $EXEC $I $R -out $O -m 4 -hitPolicy randombest -randomSeed 1 -nproc 10
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $O > $TMP1 && mv $TMP1 $O
  $ diff $O $X && echo $?
  0

Test hitpolicy leftmost
  $ O=$OUTDIR/hitpolicy_leftmost.m4
  $ X=$STDDIR/hitpolicy_leftmost.m4
  $ $EXEC $I $R -out $O -m 4 -hitPolicy leftmost -nproc 10
  [INFO]* (glob)
  [INFO]* (glob)
  $ # target is lambda x 6, leftmost -> only map to the very first x.
  $ awk '$10 > 48502 {print}' $O |wc -l 
  0
