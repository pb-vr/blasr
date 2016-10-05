Set up
  $ . $TESTDIR/setup.sh

Test input.fofn containing a new bas.h5 file. Note that the new bas.h5 file does not 
contain any /PulseData, instead contains /MultiPart/Parts.
  $ rm -f $TMP1
  $ BASFILE=/pbi/dept/secondary/siv/testdata/BlasrTestData/ctest/data/aggressiveIntervalCut/m130812_185809_42141_c100533960310000001823079711101380_s1_p0.bas.h5
  $ REFFA=/pbi/dept/secondary/siv/testdata/BlasrTestData/ctest/data/references/Ecoli_BL21_O26/sequence/Ecoli_BL21_O26.fasta
  $ $EXEC $BASFILE $REFFA --holeNumbers 1--100 --out $TMP1
  [INFO] * [blasr] started. (glob)
  [INFO] * [blasr] ended. (glob)
  $ echo $?
  0
