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

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "fmt/format.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

char *GZIP_BINARY = nullptr;

using ::testing::Eq;

class DumpAtomGZTest : public MeltTest {
    std::string dump_style = "atom";

public:
    void enable_triclinic()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command("change_box all triclinic");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void generate_dump(std::string dump_file, std::string dump_modify_options, int ntimesteps)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id all {} 1 {}", dump_style, dump_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {}", ntimesteps));
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void generate_text_and_compressed_dump(std::string text_file, std::string compressed_file,
                                           std::string compression_style,
                                           std::string dump_modify_options, int ntimesteps)
    {
        generate_text_and_compressed_dump(text_file, compressed_file, compression_style,
                                          dump_modify_options, dump_modify_options, ntimesteps);
    }

    void generate_text_and_compressed_dump(std::string text_file, std::string compressed_file,
                                           std::string compression_style,
                                           std::string text_options, std::string compressed_options, int ntimesteps)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id0 all {} 1 {}", dump_style, text_file));
        command(fmt::format("dump id1 all {} 1 {}", compression_style, compressed_file));

        if (!text_options.empty()) {
            command(fmt::format("dump_modify id0 {}", text_options));
        }

        if (!compressed_options.empty()) {
            command(fmt::format("dump_modify id1 {}", compressed_options));
        }

        command(fmt::format("run {}", ntimesteps));
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    std::string convert_compressed_to_text(std::string compressed_file)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        std::string converted_file = compressed_file.substr(0, compressed_file.find_last_of('.'));
        std::string cmdline =
            fmt::format("{} -d -c {} > {}", GZIP_BINARY, compressed_file, converted_file);
        system(cmdline.c_str());
        if (!verbose) ::testing::internal::GetCapturedStdout();
        return converted_file;
    }
};

//-------------------------------------------------------------------------------------------------
// GZ compressed files
//-------------------------------------------------------------------------------------------------

TEST_F(DumpAtomGZTest, compressed_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_run0.melt";
    auto compressed_file = "dump_gz_compressed_run0.melt.gz";

    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "", 0);

    TearDown();

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(compressed_file);

    auto converted_file = convert_compressed_to_text(compressed_file);

    ASSERT_THAT(converted_file, Eq("dump_gz_compressed_run0.melt"));
    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(compressed_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomGZTest, compressed_multi_file_run1)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_multi_file_run1_*.melt";
    auto compressed_file = "dump_gz_compressed_multi_file_run1_*.melt.gz";
    auto text_file_0     = "dump_gz_text_multi_file_run1_0.melt";
    auto text_file_1     = "dump_gz_text_multi_file_run1_1.melt";
    auto compressed_file_0 = "dump_gz_compressed_multi_file_run1_0.melt.gz";
    auto compressed_file_1 = "dump_gz_compressed_multi_file_run1_1.melt.gz";

    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "", 1);

    TearDown();

    ASSERT_FILE_EXISTS(text_file_0);
    ASSERT_FILE_EXISTS(text_file_1);
    ASSERT_FILE_EXISTS(compressed_file_0);
    ASSERT_FILE_EXISTS(compressed_file_1);

    auto converted_file_0 = convert_compressed_to_text(compressed_file_0);
    auto converted_file_1 = convert_compressed_to_text(compressed_file_1);

    ASSERT_THAT(converted_file_0, Eq("dump_gz_compressed_multi_file_run1_0.melt"));
    ASSERT_THAT(converted_file_1, Eq("dump_gz_compressed_multi_file_run1_1.melt"));
    ASSERT_FILE_EXISTS(converted_file_0);
    ASSERT_FILE_EXISTS(converted_file_1);
    ASSERT_FILE_EQUAL(text_file_0, converted_file_0);
    ASSERT_FILE_EQUAL(text_file_1, converted_file_1);

    delete_file(text_file_0);
    delete_file(text_file_1);
    delete_file(compressed_file_0);
    delete_file(compressed_file_1);
    delete_file(converted_file_0);
    delete_file(converted_file_1);
}

TEST_F(DumpAtomGZTest, compressed_multi_file_with_pad_run1)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_multi_file_pad_run1_*.melt";
    auto compressed_file = "dump_gz_compressed_multi_file_pad_run1_*.melt.gz";
    auto text_file_0     = "dump_gz_text_multi_file_pad_run1_000.melt";
    auto text_file_1     = "dump_gz_text_multi_file_pad_run1_001.melt";
    auto compressed_file_0 = "dump_gz_compressed_multi_file_pad_run1_000.melt.gz";
    auto compressed_file_1 = "dump_gz_compressed_multi_file_pad_run1_001.melt.gz";

    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "pad 3", 1);

    TearDown();

    ASSERT_FILE_EXISTS(text_file_0);
    ASSERT_FILE_EXISTS(text_file_1);
    ASSERT_FILE_EXISTS(compressed_file_0);
    ASSERT_FILE_EXISTS(compressed_file_1);

    auto converted_file_0 = convert_compressed_to_text(compressed_file_0);
    auto converted_file_1 = convert_compressed_to_text(compressed_file_1);

    ASSERT_THAT(converted_file_0, Eq("dump_gz_compressed_multi_file_pad_run1_000.melt"));
    ASSERT_THAT(converted_file_1, Eq("dump_gz_compressed_multi_file_pad_run1_001.melt"));
    ASSERT_FILE_EXISTS(converted_file_0);
    ASSERT_FILE_EXISTS(converted_file_1);
    ASSERT_FILE_EQUAL(text_file_0, converted_file_0);
    ASSERT_FILE_EQUAL(text_file_1, converted_file_1);

    delete_file(text_file_0);
    delete_file(text_file_1);
    delete_file(compressed_file_0);
    delete_file(compressed_file_1);
    delete_file(converted_file_0);
    delete_file(converted_file_1);
}

TEST_F(DumpAtomGZTest, compressed_multi_file_with_maxfiles_run1)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_multi_file_run1_*.melt";
    auto compressed_file = "dump_gz_compressed_multi_file_run1_*.melt.gz";
    auto text_file_0     = "dump_gz_text_multi_file_run1_0.melt";
    auto text_file_1     = "dump_gz_text_multi_file_run1_1.melt";
    auto text_file_2     = "dump_gz_text_multi_file_run1_2.melt";
    auto compressed_file_0 = "dump_gz_compressed_multi_file_run1_0.melt.gz";
    auto compressed_file_1 = "dump_gz_compressed_multi_file_run1_1.melt.gz";
    auto compressed_file_2 = "dump_gz_compressed_multi_file_run1_2.melt.gz";

    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "maxfiles 2", 2);

    TearDown();

    ASSERT_FILE_NOT_EXISTS(text_file_0);
    ASSERT_FILE_EXISTS(text_file_1);
    ASSERT_FILE_EXISTS(text_file_2);
    ASSERT_FILE_NOT_EXISTS(compressed_file_0);
    ASSERT_FILE_EXISTS(compressed_file_1);
    ASSERT_FILE_EXISTS(compressed_file_2);

    auto converted_file_1 = convert_compressed_to_text(compressed_file_1);
    auto converted_file_2 = convert_compressed_to_text(compressed_file_2);

    ASSERT_THAT(converted_file_1, Eq("dump_gz_compressed_multi_file_run1_1.melt"));
    ASSERT_THAT(converted_file_2, Eq("dump_gz_compressed_multi_file_run1_2.melt"));
    ASSERT_FILE_EXISTS(converted_file_1);
    ASSERT_FILE_EXISTS(converted_file_2);
    ASSERT_FILE_EQUAL(text_file_1, converted_file_1);
    ASSERT_FILE_EQUAL(text_file_2, converted_file_2);

    delete_file(text_file_1);
    delete_file(text_file_2);
    delete_file(compressed_file_1);
    delete_file(compressed_file_2);
    delete_file(converted_file_1);
    delete_file(converted_file_2);
}

TEST_F(DumpAtomGZTest, compressed_with_units_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_with_units_run0.melt";
    auto compressed_file = "dump_gz_compressed_with_units_run0.melt.gz";

    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "scale no units yes",
                                      0);

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

TEST_F(DumpAtomGZTest, compressed_with_time_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_with_time_run0.melt";
    auto compressed_file = "dump_gz_compressed_with_time_run0.melt.gz";

    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "scale no time yes",
                                      0);

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

TEST_F(DumpAtomGZTest, compressed_triclinic_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_tri_run0.melt";
    auto compressed_file = "dump_gz_compressed_tri_run0.melt.gz";

    enable_triclinic();
    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "", 0);

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

TEST_F(DumpAtomGZTest, compressed_triclinic_with_units_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_tri_with_units_run0.melt";
    auto compressed_file = "dump_gz_compressed_tri_with_units_run0.melt.gz";

    enable_triclinic();
    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "scale no units yes",
                                      0);

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

TEST_F(DumpAtomGZTest, compressed_triclinic_with_time_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_tri_with_time_run0.melt";
    auto compressed_file = "dump_gz_compressed_tri_with_time_run0.melt.gz";

    enable_triclinic();
    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "scale no time yes",
                                      0);

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

TEST_F(DumpAtomGZTest, compressed_triclinic_with_image_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_tri_with_image_run0.melt";
    auto compressed_file = "dump_gz_compressed_tri_with_image_run0.melt.gz";

    enable_triclinic();
    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "image yes", 0);

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

TEST_F(DumpAtomGZTest, compressed_modify_bad_param)
{
    command("dump id1 all atom/gz 1 dump_gz_text_modify_bad_param_run0_*.melt.gz");

    TEST_FAILURE(".*ERROR: Illegal dump_modify command: compression level must in the range of.*",
        command("dump_modify id1 compression_level 12");
    );
}

TEST_F(DumpAtomGZTest, compressed_modify_multi_bad_param)
{
    command("dump id1 all atom/gz 1 dump_gz_text_modify_bad_param_run0_*.melt.gz");

    TEST_FAILURE(".*ERROR: Illegal dump_modify command: compression level must in the range of.*",
        command("dump_modify id1 pad 3 compression_level 12");
    );
}

TEST_F(DumpAtomGZTest, compressed_modify_clevel_run0)
{
    if (!GZIP_BINARY) GTEST_SKIP();

    auto text_file       = "dump_gz_text_modify_clevel_run0.melt";
    auto compressed_file = "dump_gz_compressed_modify_clevel_run0.melt.gz";

    generate_text_and_compressed_dump(text_file, compressed_file, "atom/gz", "", "compression_level 3", 0);

    TearDown();

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(compressed_file);

    auto converted_file = convert_compressed_to_text(compressed_file);

    ASSERT_THAT(converted_file, Eq("dump_gz_compressed_modify_clevel_run0.melt"));
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

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    GZIP_BINARY = getenv("GZIP_BINARY");

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
