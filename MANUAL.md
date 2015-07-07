##NEWS:Refactored blasr has been online! Try it out!##
  
  
##Installation##

###1. Download source code###
  In order to download blasr source code from github, 
  
  + > git clone git@github.com:PacificBiosciences/blasr.git --recursive
  
  To sync your local repo with latest from github blasr master branch
  
  + > git pull -u origin master && git submodule update --init --recursive 

###2. Requirements###

  To build BLASR, you must have hdf 1.8.12 or above installed and
  configured with c++ support (you should have the library
  libhdf5_cpp.a).  If you are intalling the entire PacBio secondary
  analysis software suite, appropriate hdf libraries are already
  distributed and no configuration is necessary.  Otherwise, it is
  necessary to point two environment variables:

  + **HDF5_INC**, which points to locations of the HDF5 headers
  (e.g., hdf5.h)

  + **HDF5_LIB**, which points to the HDF5 library (e.g., hdf5*.a,
  and hdf5*.so)
  
  You may edit file blasr_git_common.mk in order sepcify HDF5_INC and HDF5_LIB.
  Alternatively, you may export them from command line.
  
  + > export HDF5_INC=path_to_your_hdf5_include && export HDF5_LIB=path_to_your_hdf5_lib

###3. Build blasr only###

  + > make blasr
  + > ./blasr

###4. Build the source tree, including other executables###

  + > make
  + > ls utils  # Frequently used executables will be under utils.


##Running BLASR##

  Typing 'blasr -h' or 'blasr -help' on the command line will give you a
  list of options.  At the least, provide a fasta, fastq, or bas.h5 file,
  and a genome.

  *Some typical use cases:*

  Align reads from reads.bas.h5 to ecoli_K12 genome, and output in SAM format.

  + > blasr reads.bas.h5  ecoli_K12.fasta -sam

  Same as above, but with soft clipping

  + > blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft

  Use multiple threads

  + > blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft -out alignments.sam -nproc 16

  Include a larger minimal match, for faster but less sensitive alignments

  + > blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft -minMatch 15

  Produce alignments in a pairwise human readable format

  + > blasr reads.bas.h5  ecoli_K12.fasta -m 0

  Use a precomputed suffix array for faster startup

  + > sawriter hg19.fasta.sa hg19.fasta #First precompute the suffix array
  + > blasr reads.bas.h5 hg19.fasta -sa hg19.fasta.sa

  Use a precomputed BWT-FM index for smaller runtime memory footprint, but slower alignments.

  + > sa2bwt hg19.fasta hg19.fasta.sa hg19.fasta.bwt
  + > blasr reads.bas.h5 hg19.fasta -bwt hg19.fasta.bwt
  

