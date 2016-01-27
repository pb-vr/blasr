
# test case headers
set( Bax2BamTest_H

    ${Bax2Bam_TestsDir}/src/TestData.h
    ${Bax2Bam_TestsDir}/src/TestUtils.h
)

# test case sources
set( Bax2BamTest_CPP

    ${Bax2Bam_TestsDir}/src/test_ccs.cpp
    ${Bax2Bam_TestsDir}/src/test_common.cpp
    ${Bax2Bam_TestsDir}/src/test_hqregions.cpp
    ${Bax2Bam_TestsDir}/src/test_polymerase.cpp
    ${Bax2Bam_TestsDir}/src/test_subreads.cpp
    ${Bax2Bam_TestsDir}/src/TestUtils.cpp
)

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

