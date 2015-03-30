Set up
  $ . $TESTDIR/setup.sh

Test blasr on ecoli.
Test blasr with -sam

# The following job takes a very long time to finish, let us use a subset of reads instead
#See $STDOUT/ecoli_v1.4.sam for 1.4 output.
# $STDOUT/ecoli_2014_03_28.sam for bug before mapQV for affineAlign/align without QV is fixed.

#  $ rm -rf $OUTDIR/ecoli.sam
#  $ $EXEC $DATDIR/ecoli.fasta $DATDIR/ecoli_reference.fasta -sam -out $OUTDIR/ecoli.sam -nproc 15
#  [INFO]* (glob)
#  [INFO]* (glob)
#
#  $ sed -n '5,$ p' $OUTDIR/ecoli.sam | sort | cut -f 1-11 > $TMP1
#  $ sed -n '5,$ p' $STDDIR/ecoli_2014_05_04.sam | sort | cut -f 1-11 > $TMP2
#  $ diff $TMP1 $TMP2
#  $ rm $TMP1 $TMP2

  $ rm -rf $OUTDIR/ecoli_subset.sam
  $ $EXEC $DATDIR/ecoli_subset.fasta $DATDIR/ecoli_reference.fasta -sam -out $OUTDIR/ecoli_subset.sam -nproc 15
  [INFO]* (glob)
  [INFO]* (glob)

  $ sed -n '5,$ p' $OUTDIR/ecoli_subset.sam | sort | cut -f 1-11 > $TMP1
  $ sed -n '5,$ p' $STDDIR/ecoli_subset_2015_03_28.sam | sort | cut -f 1-11 > $TMP2
#changelist 148101, 148080 updated read group id; 148100 updated TLEN
  $ diff $TMP1 $TMP2
  $ rm $TMP1 $TMP2
