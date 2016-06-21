.PHONY=all

SRCDIR:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
include ${CURDIR}/../../defines.mk
include ${SRCDIR}/../../rules.mk

all: ${CURDIR}/src/* ${CURDIR}/tests/src/*
	@mkdir -p ${CURDIR}/build && \
	 cd ${CURDIR}/build && \
		cmake -DBOOST_ROOT=${BOOST_ROOT} \
          -DPacBioBAM_INCLUDE_DIRS=${PBBAM_INC} \
          -DHTSLIB_INCLUDE_DIRS=${HTSLIB_INC} \
          -DPacBioBAM_LIBRARIES=${PBBAM_LIB}/libpbbam${SH_LIB_EXT} \
          -DHTSLIB_LIBRARIES=${HTSLIB_LIB}/libhts${SH_LIB_EXT} \
          -DPBDATA_INCLUDE_DIRS=${LIBPBDATA_INC} \
          -DPBDATA_LIBRARIES=${LIBPBDATA_LIB}/libpbdata${SH_LIB_EXT} \
          -DPBIHDF_INCLUDE_DIRS=${LIBPBIHDF_INC} \
          -DPBIHDF_LIBRARIES=${LIBPBIHDF_LIB}/libpbihdf${SH_LIB_EXT} \
          -DBLASR_INCLUDE_DIRS=${LIBBLASR_INC}/ \
          -DBLASR_LIBRARIES=${LIBBLASR_LIB}/libblasr${SH_LIB_EXT} \
          -DHDF5_INCLUDE_DIRS=${HDF5_INC} \
          -DHDF5_CPP_LIBRARIES=${HDF5_LIB}/libhdf5_cpp${SH_LIB_EXT} \
          -DHDF5_LIBRARIES=${HDF5_LIB}/libhdf5${SH_LIB_EXT} \
          -DBax2Bam_EXE_LINKER_FLAGS="-Wl,--no-as-needed -ldl -pthread -lrt " \
          ../ && \
		make

clean:
	@rm -rf ${CURDIR}/bin/
	@rm -rf ${CURDIR}/build
