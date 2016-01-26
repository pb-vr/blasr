# test case headers
set(Bam2BaxTest_H
    ${Bam2Bax_TestsDir}/src/TestData.h
    ${Bam2Bax_TestsDir}/src/TestConstants.h
    ${Bam2Bax_TestsDir}/src/TestUtils.h
)

set(Bam2BaxTest_CPP
    ${Bam2Bax_TestsDir}/src/test.cpp
    ${Bam2Bax_TestsDir}/src/TestUtils.cpp
    ${Bam2Bax_TestsDir}/src/test_HDFZMWWriter.cpp
    ${Bam2Bax_TestsDir}/src/test_HDFZMWWriter.cpp
    ${Bam2Bax_TestsDir}/src/test_HDFZMWMetricsWriter.cpp
    ${Bam2Bax_TestsDir}/src/test_HDFScanDataWriter.cpp
    ${Bam2Bax_TestsDir}/src/test_HDFBaseCallsWriter.cpp
    ${Bam2Bax_TestsDir}/src/test_HDFBaxWriter.cpp
    ${Bam2Bax_TestsDir}/src/test_Bam2BaxConverter.cpp
)

#set(X
#    /home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/staging/PostPrimary/bam2bax/build/src/CMakeFiles/bam2bax.dir
#)
#
#set(Bam2BaxTest_O
#    ${X}/HDFZmwWriter.cpp.o
#    ${X}/HDFRegionsWriter.cpp.o
#    ${X}/HDFScanDataWriter.cpp.o
#    ${X}/HDFBaseCallsWriter.cpp.o
#    ${X}/HDFBaxWriter.cpp.o
#)

# GoogleTest headers
set( GTest_H

    ${GTest_IncludeDir}/gtest/gtest-death-test.h
    ${GTest_IncludeDir}/gtest/gtest-message.h
    ${GTest_IncludeDir}/gtest/gtest-param-test.h
    ${GTest_IncludeDir}/gtest/gtest-printers.h
    ${GTest_IncludeDir}/gtest/gtest-spi.h
    ${GTest_IncludeDir}/gtest/gtest-test-part.h
    ${GTest_IncludeDir}/gtest/gtest-typed-test.h
    ${GTest_IncludeDir}/gtest/gtest.h
    ${GTest_IncludeDir}/gtest/gtest_pred_impl.h
    ${GTest_IncludeDir}/gtest/gtest_prod.h
    ${GTest_IncludeDir}/gtest/internal/gtest-death-test-internal.h
    ${GTest_IncludeDir}/gtest/internal/gtest-filepath.h
    ${GTest_IncludeDir}/gtest/internal/gtest-internal.h
    ${GTest_IncludeDir}/gtest/internal/gtest-linked_ptr.h
    ${GTest_IncludeDir}/gtest/internal/gtest-param-util-generated.h
    ${GTest_IncludeDir}/gtest/internal/gtest-param-util.h
    ${GTest_IncludeDir}/gtest/internal/gtest-port.h
    ${GTest_IncludeDir}/gtest/internal/gtest-string.h
    ${GTest_IncludeDir}/gtest/internal/gtest-tuple.h
    ${GTest_IncludeDir}/gtest/internal/gtest-type-util.h

    ${GTest_SourceDir}/gtest-internal-inl.h
)

# GoogleTest sources
set( GTest_CPP
    ${GTest_SourceDir}/gtest-all.cc
    ${GTest_SourceDir}/gtest_main.cc
)

