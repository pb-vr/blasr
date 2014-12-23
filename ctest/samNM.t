Set up
  $ . $TESTDIR/setup.sh

Test Sam out nm tag
  $ rm -rf $OUTDIR/read.sam
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -sam -out $OUTDIR/read.sam
  [INFO]* (glob)
  [INFO]* (glob)
  $ tail -n+5 $OUTDIR/read.sam |cut -f 21 
  NM:i:2
  NM:i:3
  NM:i:2
  NM:i:4
