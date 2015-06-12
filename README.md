#NEWS: refactored blasr has been put online. Please try it out!#

Blasr installation and maual: 

    $ more Manual.md

To pull this project from git hub to your local system:

    $ git clone git@github.com:PacificBiosciences/blasr.git blasr --recursive

To sync your code with the latest git code base:

    $ make pullfromgit -f blasr_gitp4.mk

To specify HDF5 headers and lib on your system, 

    $ edit blasr_git_common.mk

    or

    $ export HDF5_INC=path_to_your_hdf5_include && export HDF5_LIB=path_to_your_hdf5_lib

To make 'blasr' only:

    $ make blasr

To compile all tools, including blasr, pls2fasta, loadPusles, sawriter:

    $ make 

To clean all compiled tools and lib:

    $ make cleanall

To clean compiled tools without cleaning lib:

    $ make clean
