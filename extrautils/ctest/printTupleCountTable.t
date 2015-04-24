Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: printTupleCountTable.
  $ EXEC=$TESTDIR/../printTupleCountTable

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ $EXEC $OUTDIR/ecoli_tuple.table $DATDIR/ecoli_reference.fasta 
  $ echo $?
  0

  $ md5sum $OUTDIR/ecoli_tuple.table |cut -f 1 -d ' '
  3f1ae70fd009827d6d6e56050341b5df

