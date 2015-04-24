Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: samodify.
  $ EXEC=$TESTDIR/../samodify

Define tmporary files
  $ TMP1=$OUTDIR/$$.tmp.out
  $ TMP2=$OUTDIR/$$.tmp.stdout

Make OUTDIR
  $ mkdir -p $OUTDIR

  $ $EXEC $DATDIR/ecoli_reference.sa $DATDIR/ecoli_reference.fasta $OUTDIR/ecoli_reference_blt13.sa -blt 13 
  $ echo $?
  0

  $ md5sum $OUTDIR/ecoli_reference_blt13.sa | cut -f 1 -d ' '
  ac70eef5a6e03ae8177f27b3aeacc4c5
