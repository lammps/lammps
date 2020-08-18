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

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "fmt/format.h"
#include "utils.h"
#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"

#include <string>

using ::testing::Eq;

char * BINARY2TXT_BINARY = nullptr;

class DumpAtomTest : public MeltTest {
    std::string dump_style = "atom";
public:
    void enable_triclinic() {
        if (!verbose) ::testing::internal::CaptureStdout();
        command("change_box all triclinic");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void generate_dump(std::string dump_file, std::string dump_modify_options, int ntimesteps) {
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id all {} 1 {}", dump_style, dump_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {}", ntimesteps));
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void generate_text_and_binary_dump(std::string text_file, std::string binary_file, std::string dump_modify_options, int ntimesteps) {
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id0 all {} 1 {}", dump_style, text_file));
        command(fmt::format("dump id1 all {} 1 {}", dump_style, binary_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id0 {}", dump_modify_options));
            command(fmt::format("dump_modify id1 {}", dump_modify_options));
        }

        command(fmt::format("run {}", ntimesteps));
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    std::string convert_binary_to_text(std::string binary_file) {
        if (!verbose) ::testing::internal::CaptureStdout();
        std::string cmdline = fmt::format("{} {}", BINARY2TXT_BINARY, binary_file);
        system(cmdline.c_str());
        if (!verbose) ::testing::internal::GetCapturedStdout();
        return fmt::format("{}.txt", binary_file);
    }
};

TEST_F(DumpAtomTest, run0)
{
    auto dump_file = "dump_run0.melt";
    generate_dump(dump_file, "scale yes image no", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type xs ys zs"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    ASSERT_THAT(lines[9], Eq("1 1 0 0 0"));
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, format_line_run0)
{
    auto dump_file = "dump_format_line_run0.melt";
    generate_dump(dump_file, "format line \"%d %d %20.15g %g %g\" scale yes image no", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type xs ys zs"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    ASSERT_THAT(lines[9], Eq("1 1                    0 0 0"));
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, no_scale_run0)
{
    auto dump_file = "dump_no_scale_run0.melt";
    generate_dump(dump_file, "scale no", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type x y z"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, no_buffer_no_scale_run0)
{
    auto dump_file = "dump_no_buffer_no_scale_run0.melt";
    generate_dump(dump_file, "scale no", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type x y z"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, no_buffer_with_scale_run0)
{
    auto dump_file = "dump_no_buffer_with_scale_run0.melt";
    generate_dump(dump_file, "buffer no scale yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type xs ys zs"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, with_image_run0)
{
    auto dump_file = "dump_with_image_run0.melt";
    generate_dump(dump_file, "scale no image yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type x y z ix iy iz"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, with_units_run0)
{
    auto dump_file = "dump_with_units_run0.melt";
    generate_dump(dump_file, "scale no units yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 43);
    ASSERT_THAT(lines[0], Eq("ITEM: UNITS"));
    ASSERT_THAT(lines[1], Eq("lj"));
    ASSERT_THAT(lines[10], Eq("ITEM: ATOMS id type x y z"));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, with_time_run0)
{
    auto dump_file = "dump_with_time_run0.melt";
    generate_dump(dump_file, "scale no time yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 43);
    ASSERT_THAT(lines[0], Eq("ITEM: TIME"));
    ASSERT_THAT(lines[10], Eq("ITEM: ATOMS id type x y z"));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, with_units_run1)
{
    auto dump_file = "dump_with_units_run1.melt";
    generate_dump(dump_file, "scale no units yes", 1);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 84);
    ASSERT_THAT(lines[0], Eq("ITEM: UNITS"));
    ASSERT_THAT(lines[1], Eq("lj"));
    ASSERT_THAT(lines[10], Eq("ITEM: ATOMS id type x y z"));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, no_buffer_with_scale_and_image_run0)
{
    auto dump_file = "dump_no_buffer_with_scale_and_image_run0.melt";
    generate_dump(dump_file, "buffer no scale yes image yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type xs ys zs ix iy iz"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);
    delete_file(dump_file);
}
TEST_F(DumpAtomTest, triclinic_run0)
{
    auto dump_file = "dump_triclinic_run0.melt";
    enable_triclinic();
    generate_dump(dump_file, "", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type xs ys zs"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, triclinic_with_units_run0)
{
    auto dump_file = "dump_triclinic_with_units_run0.melt";
    enable_triclinic();
    generate_dump(dump_file, "units yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 43);
    ASSERT_THAT(lines[0], Eq("ITEM: UNITS"));
    ASSERT_THAT(lines[1], Eq("lj"));
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 3);
    ASSERT_THAT(lines[10], Eq("ITEM: ATOMS id type xs ys zs"));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, triclinic_with_time_run0)
{
    auto dump_file = "dump_triclinic_with_time_run0.melt";
    enable_triclinic();
    generate_dump(dump_file, "time yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 43);
    ASSERT_THAT(lines[0], Eq("ITEM: TIME"));
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 3);
    ASSERT_THAT(lines[10], Eq("ITEM: ATOMS id type xs ys zs"));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, triclinic_with_image_run0)
{
    auto dump_file = "dump_triclinic_with_image_run0.melt";
    enable_triclinic();
    generate_dump(dump_file, "image yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[4], Eq("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type xs ys zs ix iy iz"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);
    delete_file(dump_file);
}

//-------------------------------------------------------------------------------------------------
// binary formats
//-------------------------------------------------------------------------------------------------

TEST_F(DumpAtomTest, binary_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file = "dump_text_run0.melt";
    auto binary_file = "dump_binary_run0.melt.bin";

    generate_text_and_binary_dump(text_file, binary_file, "", 0);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomTest, binary_with_units_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file = "dump_text_with_units_run0.melt";
    auto binary_file = "dump_binary_with_units_run0.melt.bin";

    generate_text_and_binary_dump(text_file, binary_file, "scale no units yes", 0);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomTest, binary_with_time_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file = "dump_text_with_time_run0.melt";
    auto binary_file = "dump_binary_with_time_run0.melt.bin";

    generate_text_and_binary_dump(text_file, binary_file, "scale no time yes", 0);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomTest, binary_triclinic_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file = "dump_text_tri_run0.melt";
    auto binary_file = "dump_binary_tri_run0.melt.bin";

    enable_triclinic();
    generate_text_and_binary_dump(text_file, binary_file, "", 0);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomTest, binary_triclinic_with_units_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file = "dump_text_tri_with_units_run0.melt";
    auto binary_file = "dump_binary_tri_with_units_run0.melt.bin";

    enable_triclinic();
    generate_text_and_binary_dump(text_file, binary_file, "scale no units yes", 0);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomTest, binary_triclinic_with_time_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file = "dump_text_tri_with_time_run0.melt";
    auto binary_file = "dump_binary_tri_with_time_run0.melt.bin";

    enable_triclinic();
    generate_text_and_binary_dump(text_file, binary_file, "scale no time yes", 0);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomTest, binary_triclinic_with_image_run0)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file = "dump_text_tri_with_image_run0.melt";
    auto binary_file = "dump_binary_tri_with_image_run0.melt.bin";

    enable_triclinic();
    generate_text_and_binary_dump(text_file, binary_file, "image yes", 0);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpAtomTest, run1)
{
    auto dump_file = "dump_run1.melt";
    generate_dump(dump_file, "", 1);

    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 82);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, run2)
{
    auto dump_file = "dump_run2.melt";
    generate_dump(dump_file, "", 2);

    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 123);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, multi_file_run1)
{
    auto dump_file = "dump_run1_*.melt";
    generate_dump(dump_file, "", 1);

    ASSERT_FILE_EXISTS("dump_run1_0.melt");
    ASSERT_FILE_EXISTS("dump_run1_1.melt");
    ASSERT_EQ(count_lines("dump_run1_0.melt"), 41);
    ASSERT_EQ(count_lines("dump_run1_1.melt"), 41);
    delete_file("dump_run1_0.melt");
    delete_file("dump_run1_1.melt");
}

TEST_F(DumpAtomTest, per_processor_file_run1)
{
    auto dump_file = "dump_run1_p%.melt";
    generate_dump(dump_file, "", 1);

    ASSERT_FILE_EXISTS("dump_run1_p0.melt");
    ASSERT_EQ(count_lines("dump_run1_p0.melt"), 82);
    delete_file("dump_run1_p0.melt");
}

TEST_F(DumpAtomTest, per_processor_multi_file_run1)
{
    auto dump_file = "dump_run1_p%_*.melt";
    generate_dump(dump_file, "", 1);

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
    auto reference = "dump_ref_run0.melt";
    auto dump_file = "write_dump_atom_run0.melt";

    if (!verbose) ::testing::internal::CaptureStdout();
    command(fmt::format("dump id all atom 1 {}", reference));
    command("dump_modify id scale no units yes");
    command("run 0");
    command("write_dump all atom write_dump_atom_run*.melt modify scale no units yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS(reference);
    ASSERT_FILE_EXISTS(dump_file);

    ASSERT_FILE_EQUAL(reference, dump_file);
    delete_file(reference);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, binary_write_dump)
{
    if(!BINARY2TXT_BINARY) GTEST_SKIP();

    auto reference = "dump_run0.melt.bin";
    auto dump_file = "write_dump_atom_run0_p0.melt.bin";

    if (!verbose) ::testing::internal::CaptureStdout();
    command(fmt::format("dump id all atom 1 {}", reference));
    command("dump_modify id scale no units yes");
    command("run 0");
    command("write_dump all atom write_dump_atom_run*_p%.melt.bin modify scale no units yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS(reference);
    ASSERT_FILE_EXISTS(dump_file);

    auto reference_txt = convert_binary_to_text(reference);
    auto dump_file_txt = convert_binary_to_text(dump_file);

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
