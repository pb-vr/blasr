Set up
  $ . $TESTDIR/setup.sh

Test blasr with input bam which has:
(1) insertionQV, deletionQV, deletionTag, substitutionQV, substitutionTag
(2) insertionQV, deletionQV, deletionTag
(3) no QV
and then check if output is determined.

(1)
  $ name=iq--dq--sub
  $ infile=$DATDIR/test_bam/$name.subreads.bam
  $ outfile=$OUTDIR/$name.m4
  $ stdfile=$STDDIR/$name.m4
  $ rm -f $outfile
  $ $EXEC $infile  $DATDIR/lambda_ref.fasta -m 4 --out $outfile && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0
  $ sort $outfile > $outfile.tmp && mv $outfile.tmp $outfile
  $ diff $outfile $stdfile

(2)
  $ name=iq--dq
  $ infile=$DATDIR/test_bam/$name.subreads.bam
  $ outfile=$OUTDIR/$name.m4
  $ stdfile=$STDDIR/$name.m4
  $ rm -f $outfile
  $ $EXEC $infile  $DATDIR/lambda_ref.fasta -m 4 --out $outfile && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0
  $ sort $outfile > $outfile.tmp && mv $outfile.tmp $outfile
  $ diff $outfile $stdfile

(3)
  $ name=no--iq--dq
  $ infile=$DATDIR/test_bam/$name.subreads.bam
  $ outfile=$OUTDIR/$name.m4
  $ stdfile=$STDDIR/$name.m4
  $ rm -f $outfile
  $ $EXEC $infile  $DATDIR/lambda_ref.fasta -m 4 --out $outfile && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0
  $ sort $outfile > $outfile.tmp && mv $outfile.tmp $outfile
  $ diff $outfile $stdfile
