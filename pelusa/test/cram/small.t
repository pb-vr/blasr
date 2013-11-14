  $ export SMALL_FASTA=$TESTDIR/../data/small.fasta
  $ export OUTPUT=$TESTDIR/small.hits
  $ $TESTDIR/../../build/pelusa -q $SMALL_FASTA -t $SMALL_FASTA -k 8 > $OUTPUT
  Initializing bloom filters
  Populating bloom filters
  Pelusa collision count 0
  Querying bloom filters
  Creating worker 0
  Starting worker 0
  Finishing worker 0
  Finished pelusa
  $ cat $OUTPUT
  m130804_010516_42211_c100518912550000001823081209281361_s1_p0/15/0_8332\t\t0 (esc)
  m130804_010516_42211_c100518912550000001823081209281361_s1_p0/15/0_8332\tm130804_010516_42211_c100518912550000001823081209281361_s1_p0/15/0_8332\t6923 (esc)

