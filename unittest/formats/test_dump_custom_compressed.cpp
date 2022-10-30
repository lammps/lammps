/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

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
    DumpCustomCompressTest() : CompressedDumpTest("custom") {}
};

TEST_F(DumpCustomCompressTest, compressed_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    auto base_name       = "custom_run1.melt";
    auto text_file       = text_dump_filename(base_name);
    auto compressed_file = compressed_dump_filename(base_name);
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    if (compression_style == "custom/zstd") {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, fields, "units yes",
                                          "units yes checksum yes", 1);
    } else {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, "units yes", 1);
    }

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

TEST_F(DumpCustomCompressTest, compressed_with_time_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    auto base_name       = "with_time_custom_run1.melt";
    auto text_file       = text_dump_filename(base_name);
    auto compressed_file = compressed_dump_filename(base_name);
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    if (compression_style == "custom/zstd") {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, fields, "time yes",
                                          "time yes checksum yes", 1);
    } else {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, "time yes", 1);
    }

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

TEST_F(DumpCustomCompressTest, compressed_no_buffer_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    auto base_name       = "no_buffer_custom_run1.melt";
    auto text_file       = text_dump_filename(base_name);
    auto compressed_file = compressed_dump_filename(base_name);
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    if (compression_style == "custom/zstd") {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, fields, "buffer no",
                                          "buffer no checksum yes", 1);
    } else {
        generate_text_and_compressed_dump(text_file, compressed_file, fields, "buffer no", 1);
    }

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
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

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

TEST_F(DumpCustomCompressTest, compressed_multi_file_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    auto base_name         = "multi_file_run1_*.melt.custom";
    auto base_name_0       = "multi_file_run1_0.melt.custom";
    auto base_name_1       = "multi_file_run1_1.melt.custom";
    auto text_file         = text_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto text_file_1       = text_dump_filename(base_name_1);
    auto compressed_file   = compressed_dump_filename(base_name);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto compressed_file_1 = compressed_dump_filename(base_name_1);
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    if (compression_style == "custom/zstd") {
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

TEST_F(DumpCustomCompressTest, compressed_multi_file_with_pad_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    auto base_name         = "multi_file_pad_run1_*.melt.custom";
    auto base_name_0       = "multi_file_pad_run1_000.melt.custom";
    auto base_name_1       = "multi_file_pad_run1_001.melt.custom";
    auto text_file         = text_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto text_file_1       = text_dump_filename(base_name_1);
    auto compressed_file   = compressed_dump_filename(base_name);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto compressed_file_1 = compressed_dump_filename(base_name_1);
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

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

TEST_F(DumpCustomCompressTest, compressed_multi_file_with_maxfiles_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    auto base_name         = "multi_file_maxfiles_run1_*.melt.custom";
    auto base_name_0       = "multi_file_maxfiles_run1_0.melt.custom";
    auto base_name_1       = "multi_file_maxfiles_run1_1.melt.custom";
    auto base_name_2       = "multi_file_maxfiles_run1_2.melt.custom";
    auto text_file         = text_dump_filename(base_name);
    auto text_file_0       = text_dump_filename(base_name_0);
    auto text_file_1       = text_dump_filename(base_name_1);
    auto text_file_2       = text_dump_filename(base_name_2);
    auto compressed_file   = compressed_dump_filename(base_name);
    auto compressed_file_0 = compressed_dump_filename(base_name_0);
    auto compressed_file_1 = compressed_dump_filename(base_name_1);
    auto compressed_file_2 = compressed_dump_filename(base_name_2);
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

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

TEST_F(DumpCustomCompressTest, compressed_modify_bad_param)
{
    if (compression_style != "custom/gz") GTEST_SKIP();

    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";
    command(fmt::format("dump id1 all {} 1 {} {}", compression_style,
                        compressed_dump_filename("modify_bad_param_run0_*.melt.custom"), fields));

    TEST_FAILURE(
        ".*ERROR on proc 0: Illegal dump_modify command: Compression level must in the range of.*",
        command("dump_modify id1 compression_level 12"););
}

TEST_F(DumpCustomCompressTest, compressed_modify_multi_bad_param)
{
    if (compression_style != "custom/gz") GTEST_SKIP();

    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";
    command(fmt::format("dump id1 all {} 1 {} {}", compression_style,
                        compressed_dump_filename("modify_multi_bad_param_run0_*.melt.custom"),
                        fields));

    TEST_FAILURE(
        ".*ERROR on proc 0: Illegal dump_modify command: Compression level must in the range of.*",
        command("dump_modify id1 pad 3 compression_level 12"););
}

TEST_F(DumpCustomCompressTest, compressed_modify_clevel_run0)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    auto base_name       = "modify_clevel_run0.melt.custom";
    auto text_file       = text_dump_filename(base_name);
    auto compressed_file = compressed_dump_filename(base_name);

    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";
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
