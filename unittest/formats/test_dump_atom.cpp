/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "fmt/format.h"
#include "utils.h"
#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"

char * BINARY2TXT_BINARY = nullptr;

class DumpAtomTest : public MeltTest {
};

TEST_F(DumpAtomTest, run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run0.melt");
    command("dump_modify id scale yes image no");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_run0.melt");
    auto lines = read_lines("dump_run0.melt");
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type xs ys zs");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file("dump_run0.melt");
}

TEST_F(DumpAtomTest, no_scale_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_no_scale_run0.melt");
    command("dump_modify id scale no");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_no_scale_run0.melt");
    auto lines = read_lines("dump_no_scale_run0.melt");
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type x y z");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file("dump_no_scale_run0.melt");
}

TEST_F(DumpAtomTest, no_buffer_no_scale_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_no_buffer_no_scale_run0.melt");
    command("dump_modify id buffer no scale no");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_no_buffer_no_scale_run0.melt");
    auto lines = read_lines("dump_no_buffer_no_scale_run0.melt");
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type x y z");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file("dump_no_buffer_no_scale_run0.melt");
}

TEST_F(DumpAtomTest, no_buffer_with_scale_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_no_buffer_with_scale_run0.melt");
    command("dump_modify id buffer no scale yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_no_buffer_with_scale_run0.melt");
    auto lines = read_lines("dump_no_buffer_with_scale_run0.melt");
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type xs ys zs");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file("dump_no_buffer_with_scale_run0.melt");
}

TEST_F(DumpAtomTest, with_image_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_with_image_run0.melt");
    command("dump_modify id scale no image yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_with_image_run0.melt");

    auto lines = read_lines("dump_with_image_run0.melt");
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type x y z ix iy iz");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);

    delete_file("dump_with_image_run0.melt");
}

TEST_F(DumpAtomTest, with_units_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_with_units_run0.melt");
    command("dump_modify id scale no units yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_with_units_run0.melt");

    auto lines = read_lines("dump_with_units_run0.melt");
    ASSERT_EQ(lines.size(), 43);
    ASSERT_STREQ(lines[0].c_str(), "ITEM: UNITS");
    ASSERT_STREQ(lines[1].c_str(), "lj");
    ASSERT_STREQ(lines[10].c_str(), "ITEM: ATOMS id type x y z");
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);

    delete_file("dump_with_units_run0.melt");
}

TEST_F(DumpAtomTest, with_time_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_with_time_run0.melt");
    command("dump_modify id scale no time yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_with_time_run0.melt");

    auto lines = read_lines("dump_with_time_run0.melt");
    ASSERT_EQ(lines.size(), 43);
    ASSERT_STREQ(lines[0].c_str(), "ITEM: TIME");
    ASSERT_STREQ(lines[10].c_str(), "ITEM: ATOMS id type x y z");
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);

    delete_file("dump_with_time_run0.melt");
}

TEST_F(DumpAtomTest, with_units_run1)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_with_units_run1.melt");
    command("dump_modify id scale no units yes");
    command("run 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_with_units_run1.melt");

    auto lines = read_lines("dump_with_units_run1.melt");
    ASSERT_EQ(lines.size(), 84);
    ASSERT_STREQ(lines[0].c_str(), "ITEM: UNITS");
    ASSERT_STREQ(lines[1].c_str(), "lj");
    ASSERT_STREQ(lines[10].c_str(), "ITEM: ATOMS id type x y z");
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);

    delete_file("dump_with_units_run1.melt");
}

TEST_F(DumpAtomTest, binary_with_units_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id0 all atom 1 dump_text_with_units_run0.melt");
    command("dump id1 all atom 1 dump_binary_with_units_run0.melt.bin");
    command("dump_modify id0 scale no units yes");
    command("dump_modify id1 scale no units yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_text_with_units_run0.melt");
    ASSERT_FILE_EXISTS("dump_binary_with_units_run0.melt.bin");

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} dump_binary_with_units_run0.melt.bin", BINARY2TXT_BINARY);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_binary_with_units_run0.melt.bin.txt");
    ASSERT_FILE_EQUAL("dump_text_with_units_run0.melt", "dump_binary_with_units_run0.melt.bin.txt");
    delete_file("dump_text_with_units_run0.melt");
    delete_file("dump_binary_with_units_run0.melt.bin");
    delete_file("dump_binary_with_units_run0.melt.bin.txt");
}

TEST_F(DumpAtomTest, binary_triclinic_with_units_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("change_box all triclinic");
    command("dump id0 all atom 1 dump_text_tri_with_units_run0.melt");
    command("dump id1 all atom 1 dump_binary_tri_with_units_run0.melt.bin");
    command("dump_modify id0 scale no units yes");
    command("dump_modify id1 scale no units yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_text_tri_with_units_run0.melt");
    ASSERT_FILE_EXISTS("dump_binary_tri_with_units_run0.melt.bin");

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} dump_binary_tri_with_units_run0.melt.bin", BINARY2TXT_BINARY);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_binary_tri_with_units_run0.melt.bin.txt");
    ASSERT_FILE_EQUAL("dump_text_tri_with_units_run0.melt", "dump_binary_tri_with_units_run0.melt.bin.txt");
    delete_file("dump_text_tri_with_units_run0.melt");
    delete_file("dump_binary_tri_with_units_run0.melt.bin");
    delete_file("dump_binary_tri_with_units_run0.melt.bin.txt");
}

TEST_F(DumpAtomTest, binary_with_time_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id0 all atom 1 dump_text_with_time_run0.melt");
    command("dump id1 all atom 1 dump_binary_with_time_run0.melt.bin");
    command("dump_modify id0 scale no time yes");
    command("dump_modify id1 scale no time yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_text_with_time_run0.melt");
    ASSERT_FILE_EXISTS("dump_binary_with_time_run0.melt.bin");

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} dump_binary_with_time_run0.melt.bin", BINARY2TXT_BINARY);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_binary_with_time_run0.melt.bin.txt");
    ASSERT_FILE_EQUAL("dump_text_with_time_run0.melt", "dump_binary_with_time_run0.melt.bin.txt");
    delete_file("dump_text_with_time_run0.melt");
    delete_file("dump_binary_with_time_run0.melt.bin");
    delete_file("dump_binary_with_time_run0.melt.bin.txt");
}

TEST_F(DumpAtomTest, binary_triclinic_with_time_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("change_box all triclinic");
    command("dump id0 all atom 1 dump_text_tri_with_time_run0.melt");
    command("dump id1 all atom 1 dump_binary_tri_with_time_run0.melt.bin");
    command("dump_modify id0 scale no time yes");
    command("dump_modify id1 scale no time yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_text_tri_with_time_run0.melt");
    ASSERT_FILE_EXISTS("dump_binary_tri_with_time_run0.melt.bin");

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} dump_binary_tri_with_time_run0.melt.bin", BINARY2TXT_BINARY);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_binary_tri_with_time_run0.melt.bin.txt");
    ASSERT_FILE_EQUAL("dump_text_tri_with_time_run0.melt", "dump_binary_tri_with_time_run0.melt.bin.txt");
    delete_file("dump_text_tri_with_time_run0.melt");
    delete_file("dump_binary_tri_with_time_run0.melt.bin");
    delete_file("dump_binary_tri_with_time_run0.melt.bin.txt");
}

TEST_F(DumpAtomTest, no_buffer_with_scale_and_image_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_no_buffer_with_scale_and_image_run0.melt");
    command("dump_modify id buffer no scale yes image yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_no_buffer_with_scale_and_image_run0.melt");
    auto lines = read_lines("dump_no_buffer_with_scale_and_image_run0.melt");
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type xs ys zs ix iy iz");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);
    delete_file("dump_no_buffer_with_scale_and_image_run0.melt");
}

TEST_F(DumpAtomTest, triclinic_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();

    command("change_box all triclinic");
    command("dump id all atom 1 dump_triclinic_run0.melt");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_triclinic_run0.melt");

    auto lines = read_lines("dump_triclinic_run0.melt");
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS xy xz yz pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);

    ASSERT_EQ(lines.size(), 41);
    delete_file("dump_triclinic_run0.melt");
}

TEST_F(DumpAtomTest, triclinic_with_units_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();

    command("change_box all triclinic");
    command("dump id all atom 1 dump_triclinic_with_units_run0.melt");
    command("dump_modify id units yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_triclinic_with_units_run0.melt");

    auto lines = read_lines("dump_triclinic_with_units_run0.melt");
    ASSERT_STREQ(lines[0].c_str(), "ITEM: UNITS");
    ASSERT_STREQ(lines[1].c_str(), "lj");
    ASSERT_STREQ(lines[6].c_str(), "ITEM: BOX BOUNDS xy xz yz pp pp pp");
    ASSERT_EQ(utils::split_words(lines[7]).size(), 3);

    ASSERT_EQ(lines.size(), 43);
    delete_file("dump_triclinic_with_units_run0.melt");
}

TEST_F(DumpAtomTest, triclinic_with_time_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();

    command("change_box all triclinic");
    command("dump id all atom 1 dump_triclinic_with_time_run0.melt");
    command("dump_modify id time yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_triclinic_with_time_run0.melt");

    auto lines = read_lines("dump_triclinic_with_time_run0.melt");
    ASSERT_STREQ(lines[0].c_str(), "ITEM: TIME");
    ASSERT_STREQ(lines[6].c_str(), "ITEM: BOX BOUNDS xy xz yz pp pp pp");
    ASSERT_EQ(utils::split_words(lines[7]).size(), 3);

    ASSERT_EQ(lines.size(), 43);
    delete_file("dump_triclinic_with_time_run0.melt");
}

TEST_F(DumpAtomTest, triclinic_with_image_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("change_box all triclinic");
    command("dump id all atom 1 dump_triclinic_with_image_run0.melt");
    command("dump_modify id image yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_triclinic_with_image_run0.melt");

    auto lines = read_lines("dump_triclinic_with_image_run0.melt");
    ASSERT_EQ(lines.size(), 41);

    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS xy xz yz pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);

    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type xs ys zs ix iy iz");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);

    delete_file("dump_triclinic_with_image_run0.melt");
}

TEST_F(DumpAtomTest, binary_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id0 all atom 1 dump_text_run0.melt");
    command("dump id1 all atom 1 dump_binary_run0.melt.bin");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_text_run0.melt");
    ASSERT_FILE_EXISTS("dump_binary_run0.melt.bin");

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} dump_binary_run0.melt.bin", BINARY2TXT_BINARY);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_binary_run0.melt.bin.txt");
    ASSERT_FILE_EQUAL("dump_text_run0.melt", "dump_binary_run0.melt.bin.txt");
    delete_file("dump_text_run0.melt");
    delete_file("dump_binary_run0.melt.bin");
    delete_file("dump_binary_run0.melt.bin.txt");
}

TEST_F(DumpAtomTest, binary_triclinic_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("change_box all triclinic");
    command("dump id0 all atom 1 dump_text_tri_run0.melt");
    command("dump id1 all atom 1 dump_binary_tri_run0.melt.bin");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_text_tri_run0.melt");
    ASSERT_FILE_EXISTS("dump_binary_tri_run0.melt.bin");

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} dump_binary_tri_run0.melt.bin", BINARY2TXT_BINARY);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_binary_tri_run0.melt.bin.txt");
    ASSERT_FILE_EQUAL("dump_text_tri_run0.melt", "dump_binary_tri_run0.melt.bin.txt");
    delete_file("dump_text_tri_run0.melt");
    delete_file("dump_binary_tri_run0.melt.bin");
    delete_file("dump_binary_tri_run0.melt.bin.txt");
}

TEST_F(DumpAtomTest, binary_triclinic_with_image_run0)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("change_box all triclinic");
    command("dump id0 all atom 1 dump_text_tri_with_image_run0.melt");
    command("dump id1 all atom 1 dump_binary_tri_with_image_run0.melt.bin");
    command("dump_modify id0 image yes");
    command("dump_modify id1 image yes");
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_text_tri_with_image_run0.melt");
    ASSERT_FILE_EXISTS("dump_binary_tri_with_image_run0.melt.bin");

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} dump_binary_tri_with_image_run0.melt.bin", BINARY2TXT_BINARY);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_binary_tri_with_image_run0.melt.bin.txt");
    ASSERT_FILE_EQUAL("dump_text_tri_with_image_run0.melt",
                      "dump_binary_tri_with_image_run0.melt.bin.txt");

    auto lines = read_lines("dump_binary_tri_with_image_run0.melt.bin.txt");
    ASSERT_EQ(lines.size(), 41);

    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS xy xz yz pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);

    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type xs ys zs ix iy iz");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);

    delete_file("dump_text_tri_with_image_run0.melt");
    delete_file("dump_binary_tri_with_image_run0.melt.bin");
    delete_file("dump_binary_tri_with_image_run0.melt.bin.txt");
}

TEST_F(DumpAtomTest, run1)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run1.melt");
    command("run 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_run1.melt");
    ASSERT_EQ(count_lines("dump_run1.melt"), 82);
    delete_file("dump_run1.melt");
}

TEST_F(DumpAtomTest, run2)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run2.melt");
    command("run 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_run2.melt");
    ASSERT_EQ(count_lines("dump_run2.melt"), 123);
    delete_file("dump_run2.melt");
}

TEST_F(DumpAtomTest, multi_file_run1)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run1_*.melt");
    command("run 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_run1_0.melt");
    ASSERT_FILE_EXISTS("dump_run1_1.melt");
    ASSERT_EQ(count_lines("dump_run1_0.melt"), 41);
    ASSERT_EQ(count_lines("dump_run1_1.melt"), 41);
    delete_file("dump_run1_0.melt");
    delete_file("dump_run1_1.melt");
}

TEST_F(DumpAtomTest, per_processor_file_run1)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run1_p%.melt");
    command("run 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_run1_p0.melt");
    ASSERT_EQ(count_lines("dump_run1_p0.melt"), 82);
    delete_file("dump_run1_p0.melt");
}

TEST_F(DumpAtomTest, per_processor_multi_file_run1)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run1_p%_*.melt");
    command("run 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS("dump_run1_p0_0.melt");
    ASSERT_FILE_EXISTS("dump_run1_p0_1.melt");
    ASSERT_EQ(count_lines("dump_run1_p0_0.melt"), 41);
    ASSERT_EQ(count_lines("dump_run1_p0_1.melt"), 41);
    delete_file("dump_run1_p0_0.melt");
    delete_file("dump_run1_p0_1.melt");
}

TEST_F(DumpAtomTest, dump_modify_scale_invalid)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump.txt");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*Illegal dump_modify command.*",
                 command("dump_modify id scale true"););
}

TEST_F(DumpAtomTest, dump_modify_image_invalid)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump.txt");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*Illegal dump_modify command.*",
                 command("dump_modify id image true"););
}

TEST_F(DumpAtomTest, dump_modify_invalid)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump.txt");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*Illegal dump_modify command.*",
                 command("dump_modify id true"););
}

TEST_F(DumpAtomTest, write_dump)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run0.melt");
    command("dump_modify id scale no units yes");
    command("run 0");
    command("write_dump all atom write_dump_atom_run*.melt modify scale no units yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto reference = "dump_run0.melt";
    auto dump_file = "write_dump_atom_run0.melt";

    ASSERT_FILE_EXISTS(reference);
    ASSERT_FILE_EXISTS(dump_file);

    ASSERT_FILE_EQUAL(reference, dump_file);
    delete_file(reference);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, binary_write_dump)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dump id all atom 1 dump_run0.melt.bin");
    command("dump_modify id scale no units yes");
    command("run 0");
    command("write_dump all atom write_dump_atom_run*_p%.melt.bin modify scale no units yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    auto reference = "dump_run0.melt.bin";
    auto dump_file = "write_dump_atom_run0_p0.melt.bin";
    auto reference_txt = "dump_run0.melt.bin.txt";
    auto dump_file_txt = "write_dump_atom_run0_p0.melt.bin.txt";

    ASSERT_FILE_EXISTS(reference);
    ASSERT_FILE_EXISTS(dump_file);

    if (!verbose) ::testing::internal::CaptureStdout();
    std::string cmdline = fmt::format("{} {}", BINARY2TXT_BINARY, reference);
    system(cmdline.c_str());
    cmdline = fmt::format("{} {}", BINARY2TXT_BINARY, dump_file);
    system(cmdline.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS(reference_txt);
    ASSERT_FILE_EXISTS(dump_file_txt);

    ASSERT_FILE_EQUAL(reference_txt, dump_file_txt);
    delete_file(reference_txt);
    delete_file(dump_file_txt);
    delete_file(reference);
    delete_file(dump_file);
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

    BINARY2TXT_BINARY = getenv("BINARY2TXT_BINARY");

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
