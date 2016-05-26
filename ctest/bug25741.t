Set up
  $ . $TESTDIR/setup.sh

bug_25741, if input bas.h5 does not contain mergeQV, blasr with --printSAMQV, --nproc>1 should not write garbage 'mq' values to output.
  $ $EXEC $DATDIR/bas_wo_mergeQV.fofn $DATDIR/lambda_ref.fasta --printSAMQV --sam --clipping subread --out $OUTDIR/out_printSAMQV.sam --nproc 12 
  [INFO]* (glob)
  [INFO]* (glob)
  $ grep 'mq' $OUTDIR/out_printSAMQV.sam |wc -l
  1
