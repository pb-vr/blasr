# NEWS: refactored blasr has been put online. Please try it out!

Blasr installation and maual: 

    more Manual.md

To pull this project from git hub to your local system:

    git clone git@github.com:PacificBiosciences/blasr.git blasr --recursive

To sync your code with the latest git code base:

    git pull -u origin master && git submodule update --init --recursive 

To specify HDF5 headers and lib on your system:

    export HDF5_INCLUDE=path_to_your_hdf5_include && export HDF5_LIB=path_to_your_hdf5_lib

To configure:

    ./configure.py --no-pbbam

or with HDF5 directories:

    ./configure.py --no-pbbam HDF5_INCLUDE=... HDF5_LIB=...

To make 'blasr' only:

    make blasr

To compile all tools, including blasr, pls2fasta, loadPulses, sawriter:

    make

To test (with **cram** installed):

    #make cramtests
    make cramfast
    ## Currently:
    ## Ran 22 tests, 0 skipped, 4 failed.

To clean all compiled tools and lib:

    make cleanall

To clean compiled tools without cleaning lib:

    make clean
