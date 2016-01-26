#!/usr/env/python
import os
import os.path as op
from pbcore.util.Process import backticks
from pbcore.io import readFofn


"""
Convert bam files to bax.h5 and test
the converted bax.h5 file through the
2.3 smrtpipe resequencing pipeline.
"""

smrtwrap = "/mnt/software/s/smrtanalysis/2.3.0.p4/smrtcmds/bin/smrtwrap"

def c_dir():
    """return path to current directory"""
    return op.dirname(op.abspath(__file__))

def execute(cmd):
    """Execute cmd and raise run time error if it failed."""
    print "CMD: %s" % cmd
    o, c, m = backticks(cmd)
    if c != 0:
        raise RuntimeError("%s failed. %s" % (cmd, m))

def mkdir(path):
    """Create output dir"""
    cmd = "rm -rf %s && mkdir %s" % (path, path)
    execute(cmd)

def bax2bam_path():
    """Return path to bax2bam"""
    cmd = "which bax2bam"
    o, c, m = backticks(cmd)
    if c != 0:
        raise RuntimeError ("could not find bax2bam")
    else:
        return o[0]

def bam2bax_path():
    """Return path to bam2bax"""
    p = os.getenv('BAM2BAX', "%s/../../bin/bam2bax" % c_dir())
    print p
    if op.exists(p):
        return p
    else:
        raise IOError ("Unable to find bam2bax %s" % p)

def sr_bam_path(prefix):
    """return path to subreads.bam given prefix"""
    return "%s.subreads.bam" % prefix

def sc_bam_path(prefix):
    """return path to scraps.bam given prefix"""
    return "%s.scraps.bam" % prefix

def bax_fn(prefix):
    """return path to bax.h5 given prefix"""
    return "%s.bax.h5" % prefix

def metadata_fn(bax_fn):
    """return path to metadata.xml, which should be
    in upper directory of bax.h5"""
    return op.join(op.dirname(op.dirname(bax_fn)),
                   "%s.metadata.xml" % op.basename(bax_fn).split(".")[0])

def settings_xml():
    """return path to settings.xml"""
    return op.join(c_dir(), "settings.xml")

def prepare_bam_inputs(i_file, o_dir):
    """Prepare subreads|scraps.bam if they do not exist."""
    bam_fn, bax_fn = parse_input_file(i_file)
    if bax_fn is not None:
        return convert_bax_to_bam(bax_fn, o_dir)
    if bam_fn is not None:
        if "subreads.bam" not in bam_fn:
            raise ValueError ("%s is not a subreads.bam file" % bam_fn)
        return bam_fn.split(".subreads.bam")[0]

def parse_input_file(i_file):
    """Parse input file, get input bam or bax.h5 file."""
    bam_fn, bax_fn = None, None
    if (i_file.endswith(".bam")):
        bam_fn = i_file
    elif (i_file.endswith(".fofn")):
        fns = [f for f in readFofn(i_file)]
        if not all([f.endswith(".bax.h5") for f in fns]) or \
            len(fns) != 1:
            raise ValueError ("%s fofn should contain exactly one bax.h5 file.")
        else:
            bax_fn = fns[0]
    elif i_file.endswith(".bax.h5"):
        bax_fn = i_file
    else:
        raise ValueError ("Unsupported file format %s" % i_file)
    return bam_fn, bax_fn


def convert_bax_to_bam(bax_fn, o_dir):
    """Convert bax.h5 to bam, return prefix of bam file."""
    movie_name = op.basename(bax_fn).split(".bax.h5")[0]
    bam_prefix = op.join(o_dir, movie_name)
    cmd = "%s %s -o %s " % (bax2bam_path(), bax_fn, bam_prefix) + \
          "--subread --pulsefeatures DeletionQV,DeletionTag,InsertionQV,IPD,PulseWidth,MergeQV,SubstitutionQV,SubstitutionTag"
    execute(cmd)
    return bam_prefix

def convert_bam_to_bax(bam_prefix, o_dir, analysis_dir):
    """Call bam2bax to convert bam files to $o_dir/$analysis_dir/{movie}.bax.h5, return prefix of bax file."""
    sr_bam = sr_bam_path(bam_prefix)
    sc_bam = sc_bam_path(bam_prefix)
    movie_name = op.basename(sr_bam).split(".subreads.bam")[0]
    bax_prefix = op.join(o_dir, analysis_dir, movie_name)

    # create output directory
    cmd = "mkdir -p %s/%s" % (o_dir, analysis_dir)
    execute(cmd)

    # call bam2bax with --metadata
    cmd = "%s %s %s -o %s --metadata " % (bam2bax_path(), sr_bam, sc_bam, bax_prefix)
    execute(cmd)

    # call ls to verify that both *.bax.h5 and metadata.xml are generated
    cmd = "ls %s %s" % (bax_fn(bax_prefix), metadata_fn(bax_fn(bax_prefix)))
    execute(cmd)
    return bax_prefix

def prepare_4_resequencing(bax_prefix, o_dir):
    """Create *.fofn and *.xml from bax.h5, return xml file"""
    # Create input.fofn
    fofn = op.join(o_dir, "input.fofn")
    #with open(fofn, 'w') as f:
    #    f.write("%s" % bax_fn(bax_prefix))
    cmd = "echo %s | xargs realpath > %s" % (bax_fn(bax_prefix), fofn)
    execute(cmd)

    # Create input.xml for smrtpipe
    xml = op.join(o_dir, "input.xml")
    cmd = "%s fofnToSmrtpipeInput.py %s > %s" % (smrtwrap, fofn, xml)
    execute(cmd)

    return xml

def run_resequencing(xml, o_dir):
    cmd = "%s " % smrtwrap + \
          "smrtpipe.py --debug -D TMP=/scratch -D SHARED_DIR=/scratch --distribute " + \
          "--params=%s " % (settings_xml()) + \
          "--output=%s " % (o_dir) + \
          "xml:%s " % (xml) + \
          "2>%s/err 1>%s/out" % (o_dir, o_dir)
    execute(cmd)

def run(i_file, o_dir):
    """
    """
    mkdir(o_dir)

    bam_prefix = prepare_bam_inputs(i_file, o_dir)

    bax_prefix = convert_bam_to_bax(bam_prefix, o_dir, analysis_dir="Analysis_Results")

    xml        = prepare_4_resequencing(bax_prefix, o_dir)

    run_resequencing(xml, o_dir)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print "Usage: %s input.subreads.bam output_dir" % (op.basename(__file__))
        print "       Convert bam files to bax.h5 and run smrtpipe resequencing."
        print "Usage: %s input.fofn output_dir" % (op.basename(__file__))
        print "       Convert input.fofn's bax.h5 files to bam, then converts bam to bax, and finally run smrtpipe resequencing."
        exit(1);

    print "current directory = %s " % c_dir()

    run(i_file=sys.argv[1], o_dir=sys.argv[2])
