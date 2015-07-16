#NEWS: refactored blasr has been put online. Please try it out!#

Blasr installation and maual: 

    $ more Manual.md

To pull this project from git hub to your local system:

    $ git clone git@github.com:PacificBiosciences/blasr.git blasr --recursive

To sync your code with the latest git code base:

    $ git pull -u origin master && git submodule update --init --recursive 

To specify HDF5 headers and lib on your system, 

    $ edit blasr_git_common.mk

    or

    $ export HDF5_INC=path_to_your_hdf5_include && export HDF5_LIB=path_to_your_hdf5_lib

To make 'blasr' only:

    $ make -f yli.makefile blasr

To compile all tools, including blasr, pls2fasta, loadPulses, sawriter:

    $ make -f yli.makefile

To clean all compiled tools and lib:

    $ make -f yli.makefile cleanall

To clean compiled tools without cleaning lib:

    $ make -f yli.makefile clean
