
Simple test to make sure bax2bam runs properly.

  $ cd $WORKSPACE
  $ cd smrtanalysis/_output/modulebuilds/bioinformatics/staging/PostPrimary/bax2bam/_output/install/binwrap-build
  $ rm -f tst1.*.bam
  $ ./bax2bam -o tst1 /pbi/dept/secondary/siv/testdata/SA3-RS/lambda/2372215/0007_tiny/Analysis_Results/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.bax.h5
