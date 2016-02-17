########################################################################
# CMake build script for Bax2Bam executable.
########################################################################

project(Bax2Bam)
cmake_minimum_required(VERSION 2.8)

# project version
set(Bax2Bam_MAJOR_VERSION 0)
set(Bax2Bam_MINOR_VERSION 0)
set(Bax2Bam_PATCH_VERSION 7)
set(Bax2Bam_VERSION
  "${Bax2Bam_MAJOR_VERSION}.${Bax2Bam_MINOR_VERSION}.${Bax2Bam_PATCH_VERSION}"
)

# build-time options
option(Bax2BAM_build_tests "Build Bax2BAM's unit tests." ON)

# main project paths
set(Bax2Bam_RootDir       ${Bax2Bam_SOURCE_DIR})
set(Bax2Bam_DocsDir       ${Bax2Bam_RootDir}/docs)
set(Bax2Bam_SourceDir     ${Bax2Bam_RootDir}/src)
set(Bax2Bam_TestsDir      ${Bax2Bam_RootDir}/tests)
set(Bax2Bam_ThirdPartyDir ${Bax2Bam_RootDir}/third-party)

if (NOT Bax2BAM_OutputDir)
    set(Bax2BAM_OutputDir ${Bax2Bam_RootDir})
endif()

set(Bax2Bam_BinDir        ${Bax2BAM_OutputDir}/bin)
file(MAKE_DIRECTORY       ${Bax2Bam_BinDir})

# shared & third-party paths
if (NOT PBDATA_ROOT_DIR)
    set(PBDATA_ROOT_DIR ${Bax2Bam_RootDir}/../../../blasr_libcpp)
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
        set(HDF5_RootDir ${Bax2Bam_RootDir}/../../../../../../prebuilt.out/hdf5/hdf5-1.8.12/ubuntu-1404)
    endif()

    set(HDF5_INCLUDE_DIRS ${HDF5_RootDir}/include)
    set(HDF5_LibDir       ${HDF5_RootDir}/lib)

    find_library(HDF5_LIBRARIES     hdf5     ${HDF5_LibDir} NO_CMAKE_SYSTEM_PATH)
    find_library(HDF5_CPP_LIBRARIES hdf5_cpp ${HDF5_LibDir} NO_CMAKE_SYSTEM_PATH)
endif()

if (NOT PacBioBAM_INCLUDE_DIRS OR NOT PacBioBAM_LIBRARIES
    OR NOT HTSLIB_INCLUDE_DIRS OR NOT HTSLIB_LIBRARIES)
    set(PacBioBAM_RootDir  ${Bax2Bam_RootDir}/../../../../../lib/cpp/pbbam)
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
set(Bax2Bam_CXX_FLAGS "-g -std=c++11 -Wall")

# quash warnings from pbdata
check_cxx_compiler_flag("-Wno-overloaded-virtual" HAS_NO_OVERLOADED_VIRTUAL)
if(HAS_NO_OVERLOADED_VIRTUAL)
    set(Bax2Bam_CXX_FLAGS "${Bax2Bam_CXX_FLAGS} -Wno-overloaded-virtual")
endif()
check_cxx_compiler_flag("-Wno-unused-private-field" HAS_NO_UNUSED_PRIVATE_FIELD)
if(HAS_NO_UNUSED_PRIVATE_FIELD)
    set(Bax2Bam_CXX_FLAGS "${Bax2Bam_CXX_FLAGS} -Wno-unused-private-field")
endif()
check_cxx_compiler_flag("-Wno-unused-variable" HAS_NO_UNUSED_VARIABLE)
if(HAS_NO_UNUSED_VARIABLE)
    set(Bax2Bam_CXX_FLAGS "${Bax2Bam_CXX_FLAGS} -Wno-unused-variable")
endif()
check_cxx_compiler_flag("-Wno-uninitialized" HAS_NO_UNINITIALIZED)
if(HAS_NO_UNINITIALIZED)
    set(Bax2Bam_CXX_FLAGS "${Bax2Bam_CXX_FLAGS} -Wno-uninitialized")
endif()

SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${BAX2BAM_LINKER_FLAGS}" )

# main exe src
add_subdirectory(src)

# testing
if(Bax2BAM_build_tests)

    enable_testing()
    
    if (NOT GTEST_SRC_DIR)
        set(GTEST_SRC_DIR ../gtest)
    endif()

    add_subdirectory(${GTEST_SRC_DIR} external/gtest/build)
    add_subdirectory(tests)
endif()

