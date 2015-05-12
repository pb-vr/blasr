Set up 
  $ . $TESTDIR/setup.sh

Set up the executable: ccsh5tobam
  $ EXEC=$TESTDIR/../ccsh5tobam
  $ SMRTWRAP=/mnt/secondary/Smrtpipe/builds/Internal_Mainline_Nightly_LastSuccessfulBuild/smrtcmds/bin/smrtwrap

  $ $SMRTWRAP python $SCRIPTDIR/test_ccsh5tobam.py $EXEC $DATDIR/test_ccsh5tobam/input.fofn $OUTDIR/test_ccsh5tobam.bam
  $ echo $?
  0
