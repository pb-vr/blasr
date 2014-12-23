Set up
  $ . $TESTDIR/setup.sh

Test -aggressiveIntervalCut.
  $ rm -f $TMP1
  $ BASFILE=/mnt/data3/vol53/2450598/0001/Analysis_Results/m130812_185809_42141_c100533960310000001823079711101380_s1_p0.bas.h5
  $ REFFA=/mnt/secondary/Smrtpipe/repository/Ecoli_BL21_O26/sequence/Ecoli_BL21_O26.fasta
  $ $EXEC $BASFILE $REFFA -holeNumbers 1-100 -out $TMP1 -aggressiveIntervalCut
  [INFO] * [blasr] started. (glob)
  [INFO] * [blasr] ended. (glob)
  $ echo $?
  0
