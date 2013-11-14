  $ export FASTA=$TESTDIR/../data/small.fasta
  $ export OUTPUT=$TESTDIR/small.hits
  $ $TESTDIR/../../build/pelusa -q $FASTA -t $FASTA -k 8 > $OUTPUT
  Initializing bloom filters
  Populating bloom filters
  Pelusa collision count 0
  Querying bloom filters
  Creating worker 0
  Starting worker 0
  Finishing worker 0
  Finished pelusa
  $ cat $OUTPUT
  15/0_8332\t15/0_8332\t6923 (esc)

