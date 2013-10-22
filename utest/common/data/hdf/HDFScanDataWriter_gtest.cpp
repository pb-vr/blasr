/*
 * ============================================================================
 *
 *       Filename:  HDFScanDataWriter_gtest.cpp
 *
 *    Description:  Test common/data/hdf/HDFScanDataWriter.h
 *
 *        Version:  1.0
 *        Created:  10/14/2013 03:04:00 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */

#include "gtest/gtest.h"
#include "data/hdf/HDFScanDataWriter.h"
#include "data/hdf/HDFFile.h"
#include "datastructures/reads/ScanData.h"

TEST(HDFScanDataWriter, Write) {
    ScanData sd;
    sd.frameRate = 100;

    HDFFile outFile;
    outFile.Open("scandata.h5", H5F_ACC_TRUNC);
    HDFScanDataWriter writer(outFile);
    writer.Write(sd);
    writer.WriteMovieName("xyz");
    writer.WriteFrameRate(70);
    writer.Close();
    outFile.Close();
}
