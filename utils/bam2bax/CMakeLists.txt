########################################################################
# CMake build script for Bam2Bax executable.
########################################################################

project(Bam2Bax CXX C)
cmake_minimum_required(VERSION 2.8)

# project version
set(Bam2Bax_MAJOR_VERSION 0)
set(Bam2Bax_MINOR_VERSION 0)
set(Bam2Bax_PATCH_VERSION 1)
set(Bam2Bax_VERSION
  "${Bam2Bax_MAJOR_VERSION}.${Bam2Bax_MINOR_VERSION}.${Bam2Bax_PATCH_VERSION}"
)

# build-time options
option(Bam2Bax_build_tests "Build Bam2Bax's unit tests." ON)

# main project paths
set(Bam2Bax_RootDir       ${Bam2Bax_SOURCE_DIR})
set(Bam2Bax_DocsDir       ${Bam2Bax_RootDir}/docs)
set(Bam2Bax_SourceDir     ${Bam2Bax_RootDir}/src)
set(Bam2Bax_TestsDir      ${Bam2Bax_RootDir}/tests)
set(Bam2Bax_ThirdPartyDir ${Bam2Bax_RootDir}/third-party)

if (NOT Bam2Bax_OutputDir)
    set(Bam2Bax_OutputDir ${Bam2Bax_RootDir})
endif()

set(Bam2Bax_BinDir        ${Bam2Bax_OutputDir}/bin)
file(MAKE_DIRECTORY       ${Bam2Bax_BinDir})

# shared & third-party paths

if (NOT PBDATA_ROOT_DIR)
    set(PBDATA_ROOT_DIR ${Bam2Bax_RootDir}/../../../blasr_libcpp)
endif()

# find (existing) libraries needed by executable and tests
if (NOT BLASR_INCLUDE_DIRS OR NOT BLASR_LIBRARIES)
    find_library(BLASR_LIBRARIES    blasr    ${PBDATA_ROOT_DIR}/alignment)
    set(BLASR_INCLUDE_DIRS ${PBDATA_ROOT_DIR}/alignment)
endif()

if (NOT PBIHDF_INCLUDE_DIRS OR NOT PBIHDF_LIBRARIES)
    find_library(PBIHDF_LIBRARIES   pbihdf   ${PBDATA_ROOT_DIR}/hdf)
    set(PBIHDF_INCLUDE_DIRS ${PBDATA_ROOT_DIR}/hdf)
endif()

if (NOT PBDATA_INCLUDE_DIRS OR NOT PBDATA_LIBRARIES)
    find_library(PBDATA_LIBRARIES   pbdata   ${PBDATA_ROOT_DIR}/pbdata)
    set(PBDATA_INCLUDE_DIRS ${PBDATA_ROOT_DIR}/pbdata)
endif()

if (NOT HDF5_INCLUDE_DIRS OR NOT HDF5_LIBRARIES)
    if (NOT HDF5_RootDir)
        set(HDF5_RootDir ${Bam2Bax_RootDir}/../../../../../../prebuilt.out/hdf5/hdf5-1.8.12/ubuntu-1404)
    endif()

    set(HDF5_INCLUDE_DIRS ${HDF5_RootDir}/include)
    set(HDF5_LibDir       ${HDF5_RootDir}/lib)

    find_library(HDF5_LIBRARIES     hdf5     ${HDF5_LibDir} NO_CMAKE_SYSTEM_PATH)
    find_library(HDF5_CPP_LIBRARIES hdf5_cpp ${HDF5_LibDir} NO_CMAKE_SYSTEM_PATH)
endif()

if (NOT PacBioBAM_INCLUDE_DIRS OR NOT PacBioBAM_LIBRARIES
    OR NOT HTSLIB_INCLUDE_DIRS OR NOT HTSLIB_LIBRARIES)
    set(PacBioBAM_RootDir ${Bam2Bax_RootDir}/../../../../../lib/cpp/pbbam)
    add_subdirectory(${PacBioBAM_RootDir} external/build/pbbam)
endif()

if (NOT Boost_INCLUDE_DIRS)
    find_package(Boost REQUIRED)
endif()

if (NOT ZLIB_LIBRARIES OR NOT ZLIB_INCLUDE_DIRS)
    find_package(ZLIB REQUIRED)
endif()

# shared CXX flags for src & tests
include(CheckCXXCompilerFlag)
set(Bam2Bax_CXX_FLAGS "-g -std=c++11 -Wall")


# quash warnings from pbdata
check_cxx_compiler_flag("-Wno-overloaded-virtual" HAS_NO_OVERLOADED_VIRTUAL)
if(HAS_NO_OVERLOADED_VIRTUAL)
    set(Bam2Bax_CXX_FLAGS "${Bam2Bax_CXX_FLAGS} -Wno-overloaded-virtual")
endif()
check_cxx_compiler_flag("-Wno-unused-private-field" HAS_NO_UNUSED_PRIVATE_FIELD)
if(HAS_NO_UNUSED_PRIVATE_FIELD)
    set(Bam2Bax_CXX_FLAGS "${Bam2Bax_CXX_FLAGS} -Wno-unused-private-field")
endif()
check_cxx_compiler_flag("-Wno-unused-variable" HAS_NO_UNUSED_VARIABLE)
if(HAS_NO_UNUSED_VARIABLE)
    set(Bam2Bax_CXX_FLAGS "${Bam2Bax_CXX_FLAGS} -Wno-unused-variable")
endif()
check_cxx_compiler_flag("-Wno-uninitialized" HAS_NO_UNINITIALIZED)
if(HAS_NO_UNINITIALIZED)
    set(Bam2Bax_CXX_FLAGS "${Bam2Bax_CXX_FLAGS} -Wno-uninitialized")
endif()
check_cxx_compiler_flag("-Wunused-but-set-variable" HAS_NO_UNUSED_BUT_SET_VARIABLE)
if(HAS_NO_UNUSED_BUT_SET_VARIABLE)
    set(Bam2Bax_CXX_FLAGS "${Bam2Bax_CXX_FLAGS} -Wunused-but-set-variable")
endif()

SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${Bam2Bax_LINKER_FLAGS}" )

# main exe src
add_subdirectory(src)

# testing
if(Bam2Bax_build_tests)

    enable_testing()
    
    if (NOT GTEST_SRC_DIR)
        set(GTEST_SRC_DIR ../gtest)
    endif()

    add_subdirectory(${GTEST_SRC_DIR} external/gtest/build)
    add_subdirectory(tests)
endif()

