Set up
  $ . $TESTDIR/setup.sh

bug_25328, unrolled resequencing test 
  $ INFA=$DATDIR/bug_25328_zmw_38131.fasta
  $ REF=$DATDIR/All4mers_circular_72x_l50256.fasta
  $ OUTFA=$OUTDIR/bug_25328.m4
  $ $EXEC $INFA $REF -bestn 1 -nCandidates 1 -forwardOnly -maxMatch 14 -m 4 -out $OUTFA
  [INFO]* (glob)
  [INFO]* (glob)

  $ awk '$7-$6 >= 15000' $OUTFA |wc -l
  1
