/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/utils.h"
#include "compressed_dump_test.h"
#include "fmt/format.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>


using ::testing::Eq;

class DumpLocalCompressTest : public CompressedDumpTest {
public:
    DumpLocalCompressTest() : CompressedDumpTest("local") {
    }
};

TEST_F(DumpLocalCompressTest, compressed_run0)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("compute comp all pair/local dist eng");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto base_name         = "run*.melt.local";
    auto base_name_0       = "run0.melt.local";
    auto text_files        = text_dump_filename(base_name);
    auto compressed_files  = compressed_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto fields            = "index c_comp[1]";

    generate_text_and_compressed_dump(text_files, compressed_files, fields, "", 0);

    TearDown();

    ASSERT_FILE_EXISTS(text_file_0);
    ASSERT_FILE_EXISTS(compressed_file_0);

    auto converted_file_0 = convert_compressed_to_text(compressed_file_0);

    ASSERT_FILE_EXISTS(converted_file_0);
    ASSERT_FILE_EQUAL(text_file_0, converted_file_0);
    delete_file(text_file_0);
    delete_file(compressed_file_0);
    delete_file(converted_file_0);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " (gz|zstd)\n\n" << std::endl;
        return 1;
    }

    if(strcmp(argv[1], "gz") == 0) {
        COMPRESS_SUFFIX = "gz";
        COMPRESS_EXTENSION = "gz";
    } else if(strcmp(argv[1], "zstd") == 0) {
        COMPRESS_SUFFIX = "zstd";
        COMPRESS_EXTENSION = "zst";
    } else {
        std::cerr << "usage: " << argv[0] << " (gz|zstd)\n\n" << std::endl;
        return 1;
    }

    COMPRESS_BINARY = getenv("COMPRESS_BINARY");

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
