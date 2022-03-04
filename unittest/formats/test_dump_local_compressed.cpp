/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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
    DumpLocalCompressTest() : CompressedDumpTest("local") {}

    void SetUp() override
    {
        CompressedDumpTest::SetUp();

        BEGIN_HIDE_OUTPUT();
        command("compute comp all pair/local dist eng");
        END_HIDE_OUTPUT();
    }
};

TEST_F(DumpLocalCompressTest, compressed_run0)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name         = "run*.melt.local";
    auto base_name_0       = "run0.melt.local";
    auto text_files        = text_dump_filename(base_name);
    auto compressed_files  = compressed_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto fields            = "index c_comp[1]";

    if (compression_style == "local/zstd") {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, fields, "",
                                          "checksum yes", 0);
    } else {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, "", 0);
    }

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

TEST_F(DumpLocalCompressTest, compressed_no_buffer_run0)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name         = "no_buffer_run*.melt.local";
    auto base_name_0       = "no_buffer_run0.melt.local";
    auto text_files        = text_dump_filename(base_name);
    auto compressed_files  = compressed_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto fields            = "index c_comp[1]";

    if (compression_style == "local/zstd") {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, fields, "buffer no",
                                          "buffer no checksum yes", 0);
    } else {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, "buffer no", 0);
    }

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

TEST_F(DumpLocalCompressTest, compressed_with_time_run0)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name         = "with_time_run*.melt.local";
    auto base_name_0       = "with_time_run0.melt.local";
    auto text_files        = text_dump_filename(base_name);
    auto compressed_files  = compressed_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto fields            = "index c_comp[1]";

    if (compression_style == "local/zstd") {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, fields, "time yes",
                                          "time yes checksum yes", 0);
    } else {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, "time yes", 0);
    }

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

TEST_F(DumpLocalCompressTest, compressed_with_units_run0)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name         = "with_units_run*.melt.local";
    auto base_name_0       = "with_units_run0.melt.local";
    auto text_files        = text_dump_filename(base_name);
    auto compressed_files  = compressed_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto fields            = "index c_comp[1]";

    if (compression_style == "local/zstd") {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, fields, "units yes",
                                          "units yes checksum yes", 0);
    } else {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, "units yes", 0);
    }

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

TEST_F(DumpLocalCompressTest, compressed_triclinic_run0)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();
    enable_triclinic();

    auto base_name         = "triclinic_run*.melt.local";
    auto base_name_0       = "triclinic_run0.melt.local";
    auto text_files        = text_dump_filename(base_name);
    auto compressed_files  = compressed_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto fields            = "index c_comp[1]";

    if (compression_style == "local/zstd") {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, fields, "",
                                          "checksum yes", 0);
    } else {
        generate_text_and_compressed_dump(text_files, compressed_files, fields, "", 0);
    }

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

TEST_F(DumpLocalCompressTest, compressed_multi_file_run1)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name         = "multi_file_run1_*.melt.local";
    auto base_name_0       = "multi_file_run1_0.melt.local";
    auto base_name_1       = "multi_file_run1_1.melt.local";
    auto text_file         = text_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto text_file_1       = text_dump_filename(base_name_1);
    auto compressed_file   = compressed_dump_filename(base_name);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto compressed_file_1 = compressed_dump_filename(base_name_1);
    auto fields            = "index c_comp[1]";

    if (compression_style == "local/zstd") {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, fields, "",
                                          "checksum no", 1);
    } else {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, "", 1);
    }

    TearDown();

    auto converted_file_0 = convert_compressed_to_text(compressed_file_0);
    auto converted_file_1 = convert_compressed_to_text(compressed_file_1);

    ASSERT_THAT(converted_file_0, Eq(converted_dump_filename(base_name_0)));
    ASSERT_THAT(converted_file_1, Eq(converted_dump_filename(base_name_1)));
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

TEST_F(DumpLocalCompressTest, compressed_multi_file_with_pad_run1)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name         = "multi_file_pad_run1_*.melt.local";
    auto base_name_0       = "multi_file_pad_run1_000.melt.local";
    auto base_name_1       = "multi_file_pad_run1_001.melt.local";
    auto text_file         = text_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto text_file_1       = text_dump_filename(base_name_1);
    auto compressed_file   = compressed_dump_filename(base_name);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto compressed_file_1 = compressed_dump_filename(base_name_1);
    auto fields            = "index c_comp[1]";

    generate_text_and_compressed_dump(text_file, compressed_file, fields, "pad 3", 1);

    TearDown();

    ASSERT_FILE_EXISTS(text_file_0);
    ASSERT_FILE_EXISTS(text_file_1);
    ASSERT_FILE_EXISTS(compressed_file_0);
    ASSERT_FILE_EXISTS(compressed_file_1);

    auto converted_file_0 = convert_compressed_to_text(compressed_file_0);
    auto converted_file_1 = convert_compressed_to_text(compressed_file_1);

    ASSERT_THAT(converted_file_0, Eq(converted_dump_filename(base_name_0)));
    ASSERT_THAT(converted_file_1, Eq(converted_dump_filename(base_name_1)));
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

TEST_F(DumpLocalCompressTest, compressed_multi_file_with_maxfiles_run1)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name         = "multi_file_maxfiles_run1_*.melt.local";
    auto base_name_0       = "multi_file_maxfiles_run1_0.melt.local";
    auto base_name_1       = "multi_file_maxfiles_run1_1.melt.local";
    auto base_name_2       = "multi_file_maxfiles_run1_2.melt.local";
    auto text_file         = text_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto text_file_1       = text_dump_filename(base_name_1);
    auto text_file_2       = text_dump_filename(base_name_2);
    auto compressed_file   = compressed_dump_filename(base_name);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto compressed_file_1 = compressed_dump_filename(base_name_1);
    auto compressed_file_2 = compressed_dump_filename(base_name_2);
    auto fields            = "index c_comp[1]";

    generate_text_and_compressed_dump(text_file, compressed_file, fields, "maxfiles 2", 2);

    TearDown();

    ASSERT_FILE_NOT_EXISTS(text_file_0);
    ASSERT_FILE_EXISTS(text_file_1);
    ASSERT_FILE_EXISTS(text_file_2);
    ASSERT_FILE_NOT_EXISTS(compressed_file_0);
    ASSERT_FILE_EXISTS(compressed_file_1);
    ASSERT_FILE_EXISTS(compressed_file_2);

    auto converted_file_1 = convert_compressed_to_text(compressed_file_1);
    auto converted_file_2 = convert_compressed_to_text(compressed_file_2);

    ASSERT_THAT(converted_file_1, Eq(converted_dump_filename(base_name_1)));
    ASSERT_THAT(converted_file_2, Eq(converted_dump_filename(base_name_2)));
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

TEST_F(DumpLocalCompressTest, compressed_modify_bad_param)
{
    if (compression_style != "local/gz") GTEST_SKIP();

    auto fields = "index c_comp[1]";

    BEGIN_HIDE_OUTPUT();
    command(fmt::format("dump id1 all {} 1 {} {}", compression_style,
                        compressed_dump_filename("modify_bad_param_run0_*.melt.local"), fields));
    END_HIDE_OUTPUT();

    TEST_FAILURE(
        ".*ERROR on proc 0: Illegal dump_modify command: Compression level must in the range of.*",
        command("dump_modify id1 compression_level 12"););
}

TEST_F(DumpLocalCompressTest, compressed_modify_multi_bad_param)
{
    if (compression_style != "local/gz") GTEST_SKIP();

    auto fields = "index c_comp[1]";

    BEGIN_HIDE_OUTPUT();
    command(fmt::format("dump id1 all {} 1 {} {}", compression_style,
                        compressed_dump_filename("modify_multi_bad_param_run0_*.melt.local"),
                        fields));
    END_HIDE_OUTPUT();

    TEST_FAILURE(
        ".*ERROR on proc 0: Illegal dump_modify command: Compression level must in the range of.*",
        command("dump_modify id1 pad 3 compression_level 12"););
}

TEST_F(DumpLocalCompressTest, compressed_modify_clevel_run0)
{
    if (!COMPRESS_BINARY) GTEST_SKIP();

    auto base_name       = "modify_clevel_run0.melt.local";
    auto text_file       = text_dump_filename(base_name);
    auto compressed_file = compressed_dump_filename(base_name);
    auto fields          = "index c_comp[1]";

    generate_text_and_compressed_dump(text_file, compressed_file, fields, fields, "",
                                      "compression_level 3", 0);

    TearDown();

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(compressed_file);

    auto converted_file = convert_compressed_to_text(compressed_file);

    ASSERT_THAT(converted_file, Eq(converted_dump_filename(base_name)));
    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(compressed_file);
    delete_file(converted_file);
}
