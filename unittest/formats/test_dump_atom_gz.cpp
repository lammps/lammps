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
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id0 all {} 1 {}", dump_style, text_file));
        command(fmt::format("dump id1 all {} 1 {}", compression_style, compressed_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id0 {}", dump_modify_options));
            command(fmt::format("dump_modify id1 {}", dump_modify_options));
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
