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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Bam2Bax_CXX_FLAGS}")

add_executable(test_bam2bax ${SOURCES})

set_target_properties(test_bam2bax PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${Bam2Bax_TestsDir}/bin
)
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
)

add_test(
    NAME UnitTests 
    WORKING_DIRECTORY ${Bam2Bax_TestsDir}/bin
    COMMAND test_bam2bax 
)
