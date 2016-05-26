Create Test Files And Run Tests
First convert *.bax.h5 to bam using bax2bam
Next convert *.bam back to *.bax.h5 using bam2bax.
Then convert generated bax.h5 to bam again using bax2bam
Finally, compare whether bam files are identical.

  $ BAX2BAM=/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/smrtcmds/bin/bax2bam
  $ BLASR=/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/smrtcmds/bin/blasr
  $ . /mnt/software/Modules/current/init/sh
  $ module load samtools
  $ SAMTOOLS=samtools

BAM2BAX=? MUST SET UP PATH TO BAM2BAX

  $ I_PREFIX=m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.1
  $ I_BAX_H5=$I_PREFIX.bax.h5
  $ I_SR_BAM=$I_PREFIX.subreads.bam
  $ I_SC_BAM=$I_PREFIX.scraps.bam
  $ I_SR_SAM=$I_PREFIX.subreads.sam
  $ I_SC_SAM=$I_PREFIX.scraps.sam

  $ O_DIR=Analysis_Results
  $ O_PREFIX=$O_DIR/${I_PREFIX}
  $ O_BAX_H5=$O_PREFIX.bax.h5
  $ O_SR_BAM=$O_PREFIX.subreads.bam
  $ O_SC_BAM=$O_PREFIX.scraps.bam
  $ O_SR_SAM=$O_PREFIX.subreads.sam
  $ O_SC_SAM=$O_PREFIX.scraps.sam
  $ O_META_XML=`echo $I_PREFIX | cut -f 1 -d '.'`.metadata.xml

Clean
  $ rm -f *.bam *.sam *.tmp $O_PREFIX.* $I_PREFIX.* ../${I_PREFIX}.metadata.xml

Copy input bax.h5
  $ cp /pbi/dept/secondary/siv/testdata/bam2bax/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.1.bax.h5 .

Make output dir
  $ mkdir -p $O_DIR

Input *.subreads.bam, *.scraps.bam
  $ $BAX2BAM $I_BAX_H5 -o $I_PREFIX --subread --losslessframes && echo $?
  0
  $ ls $I_SR_BAM $I_SC_BAM > /dev/null || echo failed to convert input bax.h5 to input bam

Output bax.h5: convert input bam to output bax.h5
  $ $BAM2BAX $I_SR_BAM $I_SC_BAM -o $O_PREFIX --metadata 1>/dev/null 2>/dev/null && echo $?
  0

Check existance of metadata.xml
  $ ls $O_META_XML >/dev/null && echo $?
  0

Output bam
echo convert output bax.h5 to bam
  $ $BAX2BAM $O_BAX_H5 -o $O_PREFIX --subread --losslessframes && echo $?
  0

To sam
  $ $SAMTOOLS view -h $I_SR_BAM -o $I_SR_SAM && cat $I_SR_SAM | grep -v '^@' > $I_SR_SAM.tmp
  $ $SAMTOOLS view -h $I_SC_BAM -o $I_SC_SAM && cat $I_SC_SAM | grep -v '^@' > $I_SC_SAM.tmp
  $ $SAMTOOLS view -h $O_SR_BAM -o $O_SR_SAM && cat $O_SR_SAM | grep -v '^@' > $O_SR_SAM.tmp
  $ $SAMTOOLS view -h $O_SC_BAM -o $O_SC_SAM && cat $O_SC_SAM | grep -v '^@' > $O_SC_SAM.tmp

diff input with output
  $ diff $I_SR_SAM.tmp $O_SR_SAM.tmp || echo I.subreads.bam and O.subreads.bam are not identical
  $ diff $I_SC_SAM.tmp $O_SC_SAM.tmp || echo I.subreads.bam and O.subreads.bam are not identical

ZMW with no HQ region
  $ $BAM2BAX /pbi/dept/secondary/siv/testdata/bam2bax/all_lq/all_lq.subreads.bam /pbi/dept/secondary/siv/testdata/bam2bax/all_lq/all_lq.scraps.bam -o Analysis_Results/all_lq 1>/dev/null 2>/dev/null && echo $?
  0

  $ h5dump -d /PulseData/Regions Analysis_Results/all_lq.bax.h5 |grep "(0,0): 47775928, 2, 0, 0, 700"
     (0,0): 47775928, 2, 0, 0, 700,

