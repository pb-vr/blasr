Set up
  $ . $TESTDIR/setup.sh

Test -useccsall with Large genome.
  $ BASFILE=/mnt/data3/vol53/2450530/0014/Analysis_Results/m130507_052228_42161_c100519212550000001823079909281305_s1_p0.3.bax.h5
  $ REFDIR=/mnt/secondary/Smrtpipe/repository/hg19_M_sorted/sequence
  $ REFFA=$REFDIR/hg19_M_sorted.fasta
  $ REFSA=$REFDIR/hg19_M_sorted.fasta.sa
  $ OUTFILE=$OUTDIR/intflow.m4
  $ $EXEC $BASFILE $REFFA -out $OUTFILE -m 4 -sa $REFSA -holeNumbers 109020
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTFILE > $TMP1 && sort $STDDIR/intflow_2014_06_10.m4 > $TMP2 && diff $TMP1 $TMP2 && echo $?
  0
