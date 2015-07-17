include defines.mk

# Note that, for now, HDF5_LIB is a directory here, and a file in libcpp.
# I do not know why buildctl has set that up, but we follow the mobs
# adaptors for now.

pblib: libcpp/configure.py libcpp/makefile
	cd libcpp && NOPBBAM=true HDF5_LIB=${HDF5_LIB}/libhdf5.so ./configure.py
	${MAKE} -C libcpp
