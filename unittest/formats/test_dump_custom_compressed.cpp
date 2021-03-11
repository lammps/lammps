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

class DumpCustomCompressTest : public CompressedDumpTest {
public:
    DumpCustomCompressTest() : CompressedDumpTest("custom") {
    }
};

TEST_F(DumpCustomCompressTest, compressed_run1)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name       = "custom_run1.melt";
    auto text_file       = text_dump_filename(base_name);
    auto compressed_file = compressed_dump_filename(base_name);
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    generate_text_and_compressed_dump(text_file, compressed_file, fields, "units yes", 1);

    TearDown();

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(compressed_file);

    auto converted_file = convert_compressed_to_text(compressed_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(compressed_file);
    delete_file(converted_file);
}

TEST_F(DumpCustomCompressTest, compressed_triclinic_run1)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name       = "custom_tri_run1.melt";
    auto text_file       = text_dump_filename(base_name);
    auto compressed_file = compressed_dump_filename(base_name);
    auto fields          = "id type proc x y z xs ys zs xsu ysu zsu vx vy vz fx fy fz";

    enable_triclinic();

    generate_text_and_compressed_dump(text_file, compressed_file, fields, "units yes", 1);

    TearDown();

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(compressed_file);

    auto converted_file = convert_compressed_to_text(compressed_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(compressed_file);
    delete_file(converted_file);
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
