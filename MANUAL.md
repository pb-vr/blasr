##Installation##

###1. Requirements###

  To build BLASR, you must have hdf 1.8.12 or above installed and
  configured with c++ support (you should have the library
  libhdf5_cpp.a).  If you are intalling the entire PacBio secondary
  analysis software suite, appropriate hdf libraries are already
  distributed and no configuration is necessary.  Otherwise, it is
  necessary to point an environment variables, **HDF5_ROOT**, which 
  points to locations of the HDF5 libraries. ${HDF5_ROOT} must contain
  a sub-directory 'include' which contains HDF5 headers (e.g., hdf5.h)
  and a sub-directory 'lib' which contains HDF5 libraries (e.g., hdf5*.a
  and hdf5*.so).

  For example:

+    > export HDF5_ROOT=~/path_to_your_hdf5_dir_contain_include_and_lib/

###2. Build blasr only###

+    > make blasr
+    > .blasr

###3. Build the source tree, including other executables###

+    > make

###4. All the other executables will be under tools ###

+    > cd tools


##Running BLASR##

Typing blasr -h or blasr -help on the command line will give you a
list of options.  At the least, provide a fasta, fastq, or bas.h5 file,
and a genome.

*Some typical use cases:*

+    Align reads from reads.bas.h5 to ecoli_K12 genome, and output in SAM format.

    * > blasr reads.bas.h5  ecoli_K12.fasta -sam

+    Same as above, but with soft clipping

    * > blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft

+    Use multiple threads

    * > blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft -out alignments.sam -nproc 16

+    Include a larger minimal match, for faster but less sensitive alignments

    * > blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft -minMatch 15

+    Produce alignments in a pairwise human readable format

    * > blasr reads.bas.h5  ecoli_K12.fasta -m 0


+    Use a precomputed suffix array for faster startup

    * > sawriter hg19.fasta.sa hg19.fasta #First precompute the suffix array
    * > blasr reads.bas.h5 hg19.fasta -sa hg19.fasta.sa


+    Use a precomputed BWT-FM index for smaller runtime memory footprint, but slower alignments.
    * > sa2bwt hg19.fasta hg19.fasta.sa hg19.fasta.bwt
    * > blasr reads.bas.h5 hg19.fasta -bwt hg19.fasta.bwt
