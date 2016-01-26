#!/usr/bin/bash
PYTHON=/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/smrtcmds/bin/python

# Test with tiny.subreads|scraps.bam --> converted.bax.h5 --> resequencing 
$PYTHON bam2bax_resequencing.py  /pbi/dept/secondary/siv/testdata/bam2bax/tiny/m150905_153119_sherri_c100907042550000001823207504291693_s1_p0.subreads.bam tiny_output || exit 1
