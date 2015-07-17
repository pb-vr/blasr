#!/usr/bin/env bash
# This will not work within Travis until have have pre-compiled HDF5
# (or least headers?). But it shows the steps.
set -ex

./configure.py --sub --no-pbbam HDF5_INCLUDE=${HDF5_INCLUDE} HDF5_LIB=${HDF5_LIB}
make -j4 init-submodule
make --debug=b -j4 all
