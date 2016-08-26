# bax2bam

## Command-line interface

```

Usage: bax2bam [options] <input files...>

bax2bam converts the legacy PacBio basecall format (bax.h5) into the BAM
basecall format.

Options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

  Input/output files:
    movie.1.bax.h5 movie.2.bax.h5 ...
                        Input files which should be from the same movie
    --xml=STRING        DataSet XML file containing a list of movie names
    -f STRING, --fofn=STRING
                        File-of-file-names containing a list of input files
    -o STRING           Prefix of output filenames. Movie name will be used if
                        no prefix provided
    --output-xml=STRING
                        Explicit output XML name. If none provided via this arg,
                        bax2bam will use -o prefix (<prefix>.dataset.xml). If
                        that is not specified either, the output XML filename
                        will be <moviename>.dataset.xml

  Output read types (mutually exclusive):
    --subread           Output subreads (default)
    --hqregion          Output HQ regions
    --polymeraseread    Output full polymerase read
    --ccs               Output CCS sequences

  Pulse feature options:
    Configure pulse features in the output BAM. Supported features include:
        Pulse Feature:    BAM tag:  Default:
        DeletionQV        dq        Y
        DeletionTag       dt        Y
        InsertionQV       iq        Y
        IPD               ip        Y
        PulseWidth        pw        Y
        MergeQV           mq        Y
        SubstitutionQV    sq        Y
        SubstitutionTag   st        N
    If this option is used, then only those features listed will be included,
    regardless of the default state.

    --pulsefeatures=STRING
                        Comma-separated list of desired pulse features, using
                        the names listed above.
                        
    --losslessframes    Store full, 16-bit IPD/PulseWidth data, instead of
                        (default) downsampled, 8-bit encoding.

  Output BAM file type:
    --internal          Output BAMs in internal mode. Currently this indicates
                        that non-sequencing ZMWs should be included in the
                        output scraps BAM file, if applicable.

```
