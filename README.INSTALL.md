## Installation

### Download source code

* To pull this project from git hub to your local system:

        git clone git://github.com/PacificBiosciences/blasr.git blasr

* To sync your code with the latest git code base:

        git pull -u origin master && git submodule update --init

### Requirements

* To configure:

        ./configure.py --no-pbbam

* or with HDF5 directories (and note that `HDF5_LIB` is a *directory* here):

        ./configure.py --no-pbbam HDF5_INCLUDE=... HDF5_LIB=...

To build BLASR, you must have hdf 1.8.12 or above installed and
  configured with c++ support (you should have the library
  libhdf5_cpp.a).  If you are intalling the entire PacBio secondary
  analysis software suite, appropriate hdf libraries are already
  distributed and no configuration is necessary.  Otherwise, it is
  necessary to point two environment variables:

  + **HDF5_INCLUDE**, which points to directory of the HDF5 headers
  (e.g., hdf5.h)

  + **HDF5_LIB**, which points to the HDF5 library directory (e.g., hdf5*.a,
  and hdf5*.so)
  
  You may pass arguments to `configure.py` as above, or you may export them from command line:
  
        export HDF5_INC=path_to_your_hdf5_include && export HDF5_LIB=path_to_your_hdf5_lib

### Build

* To make the 'libcpp' libraries:

        make build-submodule

* To make 'blasr' only:

        make blasr

* To compile all tools, including blasr, pls2fasta, loadPulses, sawriter:

        make

  * Frequently used executables will be under utils.

* To test (with **cram** installed):

        #make cramtests
        make cramfast
        ## Currently:
        ## Ran 22 tests, 0 skipped, 4 failed.

* To clean all compiled tools and lib:

        make cleanall

* To clean compiled tools without cleaning lib:

        make clean

        make blasr
        ./blasr
