include_directories(
    .
    ${BLASR_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${PacBioBAM_INCLUDE_DIRS}
    ${PBDATA_INCLUDE_DIRS}
    ${PBDATA_ROOT_DIR}
    ${PBIHDF_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
)

set(SOURCES
    main.cpp
    Bax2Bam.h
    Bax2Bam.cpp
    CcsConverter.h
    CcsConverter.cpp
    ConverterBase.h
    HqRegionConverter.h
    HqRegionConverter.cpp
    IConverter.h
    IConverter.cpp
    OptionParser.h
    OptionParser.cpp
    PolymeraseReadConverter.h
    PolymeraseReadConverter.cpp
    Settings.h
    Settings.cpp
    SubreadConverter.h
    SubreadConverter.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Bax2Bam_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Bax2Bam_EXE_LINKER_FLAGS}")


add_executable(bax2bam ${SOURCES})
set_target_properties(bax2bam PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${Bax2Bam_BinDir}
)
target_link_libraries(bax2bam 
    ${BLASR_LIBRARIES}
    ${PBIHDF_LIBRARIES}
    ${PBDATA_LIBRARIES} 
    ${HDF5_LIBRARIES}
    ${HDF5_CPP_LIBRARIES}
    ${PacBioBAM_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
)
