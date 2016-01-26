#!/usr/bin/bash
PYTHON=/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/smrtcmds/bin/python

# Test with one bax.h5 --> subreads|scraps.bam --> converted.bax.h5 --> resequencing
$PYTHON bam2bax_resequencing.py /pbi/dept/secondary/siv/testdata/bam2bax/one_bax/input.fofn one_bax_output || exit 1
