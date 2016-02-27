file(MAKE_DIRECTORY ${Bam2Bax_TestsDir}/bin)
file(MAKE_DIRECTORY ${Bam2Bax_TestsDir}/out)

# Generate path for test data
configure_file(
    ${Bam2Bax_TestsDir}/src/TestData.h.in
    ${Bam2Bax_TestsDir}/src/TestData.h
)

include(files.cmake)
include_directories(
    ${Bam2Bax_SourceDir}
    ${BLASR_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${PacBioBAM_INCLUDE_DIRS}
    ${PBDATA_INCLUDE_DIRS}
    ${PBDATA_ROOT_DIR}
    ${PBIHDF_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${gtest_SOURCE_DIR}/include
    ${gtest_SOURCE_DIR}
)

set(SOURCES
    ${Bam2BaxTest_H}
    ${Bam2BaxTest_CPP}
)

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
check_cxx_compiler_flag("-Wunused-but-set-variable" HAS_UNUSED_BUT_SET_VARIABLE)
if(HAS_UNUSED_BUT_SET_VARIABLE)
    set(Bam2Bax_CXX_FLAGS "${Bam2Bax_CXX_FLAGS} -Wunused-but-set-variable")
endif()
# NOTE: -Wno-unused-local-typedefs used to quash clang warnings w/ Boost
check_cxx_compiler_flag("-Wno-unused-local-typedef" HAS_NO_UNUSED_LOCAL_TYPEDEF)
if(HAS_NO_UNUSED_LOCAL_TYPEDEF)
    set(Bam2Bax_CXX_FLAGS "${Bam2Bax_CXX_FLAGS} -Wno-unused-local-typedef")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Bam2Bax_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Bam2Bax_EXE_LINKER_FLAGS}")

add_executable(test_bam2bax ${SOURCES})

set_target_properties(test_bam2bax PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${Bam2Bax_TestsDir}/bin
)
if (NOT APPLE)
    set(MY_LIBRT -lrt)
else()
endif()
target_link_libraries(test_bam2bax
    gtest
    gtest_main
    ${BLASR_LIBRARIES}
    ${PBIHDF_LIBRARIES}
    ${PBDATA_LIBRARIES}
    ${HDF5_CPP_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${PacBioBAM_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT} # quirky pthreads
    ${MY_LIBRT}
)

add_test(
    NAME UnitTests 
    WORKING_DIRECTORY ${Bam2Bax_TestsDir}/bin
    COMMAND test_bam2bax 
)
