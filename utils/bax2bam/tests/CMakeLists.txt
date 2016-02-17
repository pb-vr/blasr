file(MAKE_DIRECTORY  ${Bax2Bam_TestsDir}/bin)

# Generate paths for test data
configure_file(
    ${Bax2Bam_TestsDir}/src/TestData.h.in
    ${Bax2Bam_TestsDir}/src/TestData.h
)
include(files.cmake)
include_directories(
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
#    ${Bax2BamTest_H}
)

# Grab postprimary test source files

set(SOURCES
    ${Bax2BamTest_H}
    ${Bax2BamTest_CPP}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Bax2Bam_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Bax2Bam_EXE_LINKER_FLAGS}")

add_executable(test_bax2bam ${SOURCES})
# add_executable(test_bax2bam ${Bax2Bam_TestsDir}/src/test_ccs.cpp)
set_target_properties(test_bax2bam PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${Bax2Bam_TestsDir}/bin
)
if (NOT APPLE)
    set(MY_LIBRT -lrt)
else()
endif()
target_link_libraries(test_bax2bam
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
# add_test(test_bax2bam test_bax2bam)
add_test(
    NAME UnitTests
    WORKING_DIRECTORY ${Bax2Bam_TestsDir}/bin
    COMMAND test_bax2bam
)
