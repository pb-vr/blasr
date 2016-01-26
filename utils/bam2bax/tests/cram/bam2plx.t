Create Test Files And Run Tests
Convert *.bam to *.pls.h5 using bam2plx.

Must define BAM2PLX= either as 'bam2bax --pulse' or 'bam2plx'
#  $ BAM2PLX='PATH_TO_//depot/software/smrtanalysis/bioinformatics/staging/PostPrimary/bam2bax/bin/bam2bax --pulse'
#  $ BAM2PLX='PATH_TO_//depot/software/smrtanalysis/bioinformatics/staging/PostPrimary/bam2bax/bin/bam2plx'

  $ . /mnt/software/Modules/current/init/sh
  $ module load samtools
  $ SAMTOOLS=samtools

  $ I_PREFIX=$TESTDIR/../data/tiny_bam2plx
  $ I_PL_SAM=$I_PREFIX.polymerase.sam
  $ I_PL_BAM=$I_PREFIX.polymerase.bam

  $ O_DIR=`realpath Analysis_Results`
  $ O_PREFIX=$O_DIR/`basename ${I_PREFIX}`
  $ O_PLX_H5=$O_PREFIX.plx.h5
  $ O_META_XML=$O_DIR/../`basename $I_PREFIX | cut -f 1 -d '.'`.metadata.xml

Clean
  $ rm -rf $O_PLX_H5 $O_META_XML

  $ mkdir -p $O_DIR

===================================================================
Convert polymerase bam to plx.h5
  $ $SAMTOOLS view -bS $I_PL_SAM -o $I_PL_BAM 1>/dev/null 2>/dev/null && echo $?
  0

Old bam input should be rejected
  $ $BAM2PLX $I_PL_BAM -o $O_PREFIX 2>&1 |tail -1
  ERROR:* (glob)

===================================================================
Convert subreads.bam + scraps.bam to plx.h5
  $ I_PREFIX=$TESTDIR/../data/tiny
  $ I_SR_BAM=$I_PREFIX.subreads.bam
  $ I_SC_BAM=$I_PREFIX.scraps.bam
  $ I_SR_SAM=$I_PREFIX.subreads.sam
  $ I_SC_SAM=$I_PREFIX.scraps.sam

Convert input bam to output plx.h5
  $ $SAMTOOLS view -bS $I_SR_SAM -o $I_SR_BAM 1>/dev/null 2>/dev/null && echo $?
  0
  $ $SAMTOOLS view -bS $I_SC_SAM -o $I_SC_BAM 1>/dev/null 2>/dev/null && echo $?
  0
 
Check exit status and output
  $ $BAM2PLX $I_SR_BAM $I_SC_BAM -o $O_PREFIX --metadata 1>/dev/null 2>/dev/null && echo $?
  0

Check existance of metadata.xml
  $ ls $O_META_XML >/dev/null && echo $?
  0

Check h5
  $ h5ls -r $O_PLX_H5
  /                        Group
  /PulseData               Group
  /PulseData/BaseCalls     Group
  /PulseData/BaseCalls/Basecall Dataset {15576/Inf}
  /PulseData/BaseCalls/DeletionQV Dataset {15576/Inf}
  /PulseData/BaseCalls/DeletionTag Dataset {15576/Inf}
  /PulseData/BaseCalls/InsertionQV Dataset {15576/Inf}
  /PulseData/BaseCalls/MergeQV Dataset {15576/Inf}
  /PulseData/BaseCalls/PreBaseFrames Dataset {15576/Inf}
  /PulseData/BaseCalls/PulseIndex Dataset {15576/Inf}
  /PulseData/BaseCalls/QualityValue Dataset {15576/Inf}
  /PulseData/BaseCalls/SubstitutionQV Dataset {15576/Inf}
  /PulseData/BaseCalls/SubstitutionTag Dataset {15576/Inf}
  /PulseData/BaseCalls/WidthInFrames Dataset {15576/Inf}
  /PulseData/BaseCalls/ZMW Group
  /PulseData/BaseCalls/ZMW/HoleNumber Dataset {9/Inf}
  /PulseData/BaseCalls/ZMW/HoleStatus Dataset {9/Inf}
  /PulseData/BaseCalls/ZMW/HoleXY Dataset {9/Inf, 2}
  /PulseData/BaseCalls/ZMW/NumEvent Dataset {9/Inf}
  /PulseData/BaseCalls/ZMWMetrics Group
  /PulseData/BaseCalls/ZMWMetrics/HQRegionSNR Dataset {9/Inf, 4}
  /PulseData/BaseCalls/ZMWMetrics/Productivity Dataset {9/Inf}
  /PulseData/BaseCalls/ZMWMetrics/ReadScore Dataset {9/Inf}
  /PulseData/PulseCalls    Group
  /PulseData/PulseCalls/AltLabel Dataset {15581/Inf}
  /PulseData/PulseCalls/AltLabelQV Dataset {15581/Inf}
  /PulseData/PulseCalls/Channel Dataset {15581/Inf}
  /PulseData/PulseCalls/Chi2 Dataset {15581/Inf, 4}
  /PulseData/PulseCalls/IsPulse Dataset {15581/Inf}
  /PulseData/PulseCalls/LabelQV Dataset {15581/Inf}
  /PulseData/PulseCalls/MaxSignal Dataset {15581/Inf}
  /PulseData/PulseCalls/MeanSignal Dataset {15581/Inf, 4}
  /PulseData/PulseCalls/MergeQV Dataset {15581/Inf}
  /PulseData/PulseCalls/MidSignal Dataset {15581/Inf}
  /PulseData/PulseCalls/MidStdDev Dataset {15581/Inf}
  /PulseData/PulseCalls/StartFrame Dataset {15581/Inf}
  /PulseData/PulseCalls/WidthInFrames Dataset {15581/Inf}
  /PulseData/PulseCalls/ZMW Group
  /PulseData/PulseCalls/ZMW/BaselineLevel Dataset {9/Inf, 4}
  /PulseData/PulseCalls/ZMW/BaselineSigma Dataset {9/Inf, 4}
  /PulseData/PulseCalls/ZMW/HoleNumber Dataset {9/Inf}
  /PulseData/PulseCalls/ZMW/HoleStatus Dataset {9/Inf}
  /PulseData/PulseCalls/ZMW/HoleXY Dataset {9/Inf, 2}
  /PulseData/PulseCalls/ZMW/NumEvent Dataset {9/Inf}
  /PulseData/PulseCalls/ZMW/SignalLevel Dataset {9/Inf, 4}
  /PulseData/PulseCalls/ZMW/SignalSigma Dataset {9/Inf, 4}
  /PulseData/Regions       Dataset {44/Inf, 5}
  /ScanData                Group
  /ScanData/AcqParams      Group
  /ScanData/DyeSet         Group
  /ScanData/RunInfo        Group

===================================================================
Check PulseCalls/MeanSignal
  $ h5dump -d /PulseData/PulseCalls/MeanSignal $O_PLX_H5 |sed -n '2,10p'
  DATASET "/PulseData/PulseCalls/MeanSignal" {
     DATATYPE  H5T_STD_U16LE
     DATASPACE  SIMPLE { ( 15581, 4 ) / ( H5S_UNLIMITED, 4 ) }
     DATA {
     (0,0): 0, 0, 0, 112,
     (1,0): 0, 75, 0, 0,
     (2,0): 0, 0, 0, 65,
     (3,0): 0, 60, 0, 0,
     (4,0): 0, 0, 0, 62,

Check PulseCalls/MidSignal
  $ h5dump -d /PulseData/PulseCalls/MidSignal $O_PLX_H5 |grep "(0):"
     (0): 0, 0, 0, 51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 82, 0, 82,

Check PulseCalls/WidthInFrames
  $ h5dump -d /PulseData/PulseCalls/WidthInFrames $O_PLX_H5 | grep "(0):"
     (0): 1, 2, 1, 3, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 4, 1, 3, 1,


===================================================================
Check BaseCalls attributes: SchemaRevision 
  $ h5dump -a /PulseData/BaseCalls/SchemaRevision $O_PLX_H5 |grep ":"
  *"1.1" (glob)
  $ h5dump -a /PulseData/BaseCalls/SchemaRevision $O_PLX_H5 |grep "DATATYPE"
     DATATYPE  H5T_STRING {

Check PulseCalls attributes: SchemaRevision
  $ h5dump -a /PulseData/PulseCalls/ChangeListID $O_PLX_H5 | grep ":"
  *"3.0.14.167651" (glob)
  $ h5dump -a /PulseData/PulseCalls/ChangeListID $O_PLX_H5 | grep "DATATYPE"
     DATATYPE  H5T_STRING {

$ h5dump -a /PulseData/PulseCalls/Content $O_PLX_H5 | grep ":"
* "AltLabel,AltLabelQV,Channel,Chi2,IsPulse,LabelQV,MaxSignal,MeanSignal,MergeQV,MidSignal,MidStdDev,StartFrame,WidthInFrames,uint8_t,uint8_t,uint8_t,uint16_t,uint8_t,uint8_t,uint16_t,uint16_t,uint8_t,uint16_t,uint16_t,uint32_t,uint16_t" (glob)

  $ h5dump -a /PulseData/PulseCalls/ContentStored $O_PLX_H5 | grep ":"
     (0): 9
  $ h5dump -a /PulseData/PulseCalls/ContentStored $O_PLX_H5 | grep "DATATYPE"
     DATATYPE  H5T_STD_U32LE

  $ h5dump -a /PulseData/PulseCalls/DataCreated $O_PLX_H5 | grep "DATATYPE"
     DATATYPE  H5T_STRING {

  $ h5dump -a /PulseData/PulseCalls/SchemaRevision $O_PLX_H5 | grep ":"
  *"1.1" (glob)
  $ h5dump -a /PulseData/PulseCalls/SchemaRevision $O_PLX_H5 |grep "DATATYPE"
     DATATYPE  H5T_STRING {

Check ScanData/AcqParams
  $ h5dump -a /ScanData/AcqParams/AduGain $O_PLX_H5 | grep ":"
  *1 (glob)
  $ h5dump -a /ScanData/AcqParams/AduGain $O_PLX_H5 | grep "DATATYPE"
  *H5T_IEEE_F32LE (glob)

  $ h5dump -a /ScanData/AcqParams/CameraGain $O_PLX_H5 | grep ":"
  *1 (glob)
  $ h5dump -a /ScanData/AcqParams/CameraGain $O_PLX_H5 | grep "DATATYPE"
  *H5T_IEEE_F32LE (glob)

  $ h5dump -a /ScanData/AcqParams/CameraType $O_PLX_H5 | grep ":"
  *0 (glob)
  $ h5dump -a /ScanData/AcqParams/CameraType $O_PLX_H5 | grep "DATATYPE"
  *H5T_STD_I32LE (glob)

  $ h5dump -a /ScanData/AcqParams/HotStartFrame $O_PLX_H5 | grep ":"
  * 0 (glob)
  $ h5dump -a /ScanData/AcqParams/HotStartFrame $O_PLX_H5 | grep "DATATYPE"
  *H5T_STD_U32LE (glob)

  $ h5dump -a /ScanData/AcqParams/LaserOnFrame $O_PLX_H5 | grep ":"
  * 0 (glob)
  $ h5dump -a /ScanData/AcqParams/LaserOnFrame $O_PLX_H5 | grep "DATATYPE"
  *H5T_STD_U32LE (glob)

  $ h5dump -a /ScanData/AcqParams/FrameRate $O_PLX_H5 | grep ":"
  *80.047 (glob)
  $ h5dump -a /ScanData/AcqParams/FrameRate $O_PLX_H5 | grep "DATATYPE"
  *H5T_IEEE_F32LE (glob)

Check ScanData/RunInfo
  $ h5dump -a /ScanData/RunInfo/InstrumentName $O_PLX_H5 | grep ":"
  *"sequel" (glob)
  $ h5dump -a /ScanData/RunInfo/InstrumentName $O_PLX_H5 | grep "DATATYPE"
     DATATYPE  H5T_STRING {

  $ h5dump -a /ScanData/RunInfo/PlatformId $O_PLX_H5 | grep ":"
  *4 (glob)
  $ h5dump -a /ScanData/RunInfo/PlatformId $O_PLX_H5 | grep "DATATYPE"
     DATATYPE  H5T_STD_U32LE

  $ h5dump -a /ScanData/RunInfo/PlatformName $O_PLX_H5 | grep ":"
  *"SequelAlpha" (glob)
  $ h5dump -a /ScanData/RunInfo/PlatformName $O_PLX_H5 | grep "DATATYPE"
     DATATYPE  H5T_STRING {

  $ h5dump -a /ScanData/RunInfo/MovieName $O_PLX_H5 | grep ":"
  *"tiny_bam2plx" (glob)
  $ h5dump -a /ScanData/RunInfo/MovieName $O_PLX_H5 | grep "DATATYPE"
     DATATYPE  H5T_STRING {

Check ScanData/DyeSet
  $ h5dump -a /ScanData/DyeSet/BaseMap $O_PLX_H5 |grep ":"
  *"TGCA" (glob)
  $ h5dump -a /ScanData/DyeSet/BaseMap $O_PLX_H5 |grep "DATATYPE"
     DATATYPE  H5T_STRING {

  $ h5dump -a /ScanData/DyeSet/NumAnalog $O_PLX_H5 |grep "DATATYPE"
  *H5T_STD_U16LE (glob)

===================================================================
Check exit status and output
  $ $BAM2PLX $I_SR_BAM $I_SC_BAM -o ${O_PREFIX}.atgc --baseMap ATCG 1>/dev/null 2>/dev/null && echo $?
  0

Check /ScanData/DySet/BaseMap
  $ h5dump -a /ScanData/DyeSet/BaseMap $O_PREFIX.atgc.plx.h5 |grep "ATCG" |wc -l
  1
