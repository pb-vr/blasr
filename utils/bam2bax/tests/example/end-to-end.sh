#How to Create Test Files And Run Tests#

#BAX2BAM=/home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/staging/PostPrimary/bax2bam/bin/bax2bam
#BAM2BAX=/home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/staging/PostPrimary/bam2bax/bin/bam2bax
#BLASR=/home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/ext/pi/blasr/blasr

BAX2BAM=bax2bam
BAM2BAX=../../bin/bam2bax
BLASR=blasr

I_PREFIX=m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.1
I_BAX_H5=$I_PREFIX.bax.h5
I_SR_BAM=$I_PREFIX.subreads.bam
I_SC_BAM=$I_PREFIX.scraps.bam
I_SR_SAM=$I_PREFIX.subreads.sam
I_SC_SAM=$I_PREFIX.scraps.sam

O_DIR=Analysis_Results
O_PREFIX=$O_DIR/${I_PREFIX}
O_BAX_H5=$O_PREFIX.bax.h5
O_SR_BAM=$O_PREFIX.subreads.bam
O_SC_BAM=$O_PREFIX.scraps.bam
O_SR_SAM=$O_PREFIX.subreads.sam
O_SC_SAM=$O_PREFIX.scraps.sam
O_META_XML=`echo $I_PREFIX | cut -f 1 -d '.'`.metadata.xml

# Clean
rm -f *.bam *.sam *.tmp $O_PREFIX.* $I_PREFIX.*

# Make output dir
mkdir -p $O_DIR

# Input bax.h5
#echo  cp input bax.h5 from siv
cp /pbi/dept/secondary/siv/testdata/bam2bax/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.1.bax.h5 .

# Input *.subreads.bam, *.scraps.bam
echo convert input bax.h5 to bam
$BAX2BAM $I_BAX_H5 -o $I_PREFIX --subread --losslessframes
ls $I_SR_BAM $I_SC_BAM || echo failed to convert input bax.h5 to input bam

# Output bax.h5
echo convert input bam to output bax.h5
$BAM2BAX $I_SR_BAM $I_SC_BAM -o $O_PREFIX --metadata

ls $O_BAX_H5 || echo "ERROR! $O_BAX_H5" does not exist
ls $O_META_XML || echo "ERROR! $O_META_XML" does not exist

# Output bam
echo convert output bax.h5 to bam
$BAX2BAM $O_BAX_H5 -o $O_PREFIX --subread --losslessframes

# out.subreads.sam, out.scraps.sam
samtools view -h $I_SR_BAM -o $I_SR_SAM && cat $I_SR_SAM | grep -v '^@' > $I_SR_SAM.tmp
samtools view -h $I_SC_BAM -o $I_SC_SAM && cat $I_SC_SAM | grep -v '^@' > $I_SC_SAM.tmp
samtools view -h $O_SR_BAM -o $O_SR_SAM && cat $O_SR_SAM | grep -v '^@' > $O_SR_SAM.tmp
samtools view -h $O_SC_BAM -o $O_SC_SAM && cat $O_SC_SAM | grep -v '^@' > $O_SC_SAM.tmp

echo diff input with output
diff $I_SR_SAM.tmp $O_SR_SAM.tmp || echo I.subreads.bam and O.subreads.bam are not identical
diff $I_SC_SAM.tmp $O_SC_SAM.tmp || echo I.subreads.bam and O.subreads.bam are not identical

#How to use bam2bax#
#- basic use case
#    bam2bax smrtcell.subreads.bam smrtcell.scraps.bam -o output_prefix

