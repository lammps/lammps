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

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "fmt/format.h"
#include "output.h"
#include "thermo.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

using ::testing::Eq;

char *BINARY2TXT_BINARY = nullptr;
bool verbose            = false;

class DumpAtomTest : public MeltTest {
    std::string dump_style = "atom";

public:
    void enable_triclinic()
    {
        BEGIN_HIDE_OUTPUT();
        command("change_box all triclinic");
        END_HIDE_OUTPUT();
    }

    std::string dump_filename(std::string ident)
    {
        return fmt::format("dump_{}_{}.melt", dump_style, ident);
    }

    std::string text_dump_filename(std::string ident)
    {
        return fmt::format("dump_{}_text_{}.melt", dump_style, ident);
    }

    std::string binary_dump_filename(std::string ident)
    {
        return fmt::format("dump_{}_binary_{}.melt.bin", dump_style, ident);
    }

    void generate_dump(std::string dump_file, std::string dump_modify_options, int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id all {} 1 {}", dump_style, dump_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {} post no", ntimesteps));
        END_HIDE_OUTPUT();
    }

    void continue_dump(int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("run {} pre no post no", ntimesteps));
        END_HIDE_OUTPUT();
    }

    void generate_text_and_binary_dump(std::string text_file, std::string binary_file,
                                       std::string dump_modify_options, int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id0 all {} 1 {}", dump_style, text_file));
        command(fmt::format("dump id1 all {} 1 {}", dump_style, binary_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id0 {}", dump_modify_options));
            command(fmt::format("dump_modify id1 {}", dump_modify_options));
        }

        command(fmt::format("run {} post no", ntimesteps));
        END_HIDE_OUTPUT();
    }

    std::string convert_binary_to_text(std::string binary_file)
    {
        BEGIN_HIDE_OUTPUT();
        std::string cmdline = fmt::format("{} {}", BINARY2TXT_BINARY, binary_file);
        system(cmdline.c_str());
        END_HIDE_OUTPUT();
        return fmt::format("{}.txt", binary_file);
    }
};

TEST_F(DumpAtomTest, run0)
{
    auto dump_file = dump_filename("run0");
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
    auto dump_file = dump_filename("format_line_run0");
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
    auto dump_file = dump_filename("no_scale_run0");
    generate_dump(dump_file, "scale off", 0);

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
    auto dump_file = dump_filename("no_buffer_no_scale_run0");
    generate_dump(dump_file, "buffer false scale false", 0);

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
    auto dump_file = dump_filename("no_buffer_with_scale_run0");
    generate_dump(dump_file, "buffer 0 scale 1", 0);

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
    auto dump_file = dump_filename("with_image_run0");
    generate_dump(dump_file, "scale no image on", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_THAT(lines[8], Eq("ITEM: ATOMS id type x y z ix iy iz"));
    ASSERT_EQ(utils::split_words(lines[9]).size(), 8);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, with_units_run0)
{
    auto dump_file = dump_filename("with_units_run0");
    generate_dump(dump_file, "scale false units 1", 0);

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
    auto dump_file = dump_filename("with_time_run0");
    generate_dump(dump_file, "scale off time true", 0);

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
    auto dump_file = dump_filename("with_units_run1");
    generate_dump(dump_file, "scale 0 units on", 1);

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
    auto dump_file = dump_filename("no_buffer_with_scale_and_image_run0");
    generate_dump(dump_file, "buffer 0 scale 1 image true", 0);

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
    auto dump_file = dump_filename("triclinic_run0");
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
    auto dump_file = dump_filename("triclinic_with_units_run0");
    enable_triclinic();
    generate_dump(dump_file, "units on", 0);

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
    auto dump_file = dump_filename("triclinic_with_time_run0");
    enable_triclinic();
    generate_dump(dump_file, "time on", 0);

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
    auto dump_file = dump_filename("triclinic_with_image_run0");
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
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = text_dump_filename("run0");
    auto binary_file = binary_dump_filename("run0");

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
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = text_dump_filename("with_units_run0");
    auto binary_file = binary_dump_filename("with_units_run0");

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
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = text_dump_filename("with_time_run0");
    auto binary_file = binary_dump_filename("with_time_run0");

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
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = text_dump_filename("tri_run0");
    auto binary_file = binary_dump_filename("tri_run0");

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
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = text_dump_filename("tri_with_units_run0");
    auto binary_file = binary_dump_filename("tri_with_units_run0");

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
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = text_dump_filename("tri_with_time_run0");
    auto binary_file = binary_dump_filename("tri_with_time_run0");

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
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = text_dump_filename("tri_with_image_run0");
    auto binary_file = binary_dump_filename("tri_with_image_run0");

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

TEST_F(DumpAtomTest, run1plus1)
{
    auto dump_file = dump_filename("run1plus1");
    generate_dump(dump_file, "", 1);

    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 82);
    continue_dump(1);
    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 123);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, run2)
{
    auto dump_file = dump_filename("run2");
    generate_dump(dump_file, "", 2);

    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 123);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, rerun)
{
    auto dump_file = dump_filename("rerun");
    HIDE_OUTPUT([&] {
        command("fix 1 all nve");
    });
    generate_dump(dump_file, "format line \"%d %d %20.15g %20.15g %20.15g\"", 1);
    double pe_1, pe_2, pe_rerun;
    lmp->output->thermo->evaluate_keyword("pe", &pe_1);
    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 82);
    continue_dump(1);
    lmp->output->thermo->evaluate_keyword("pe", &pe_2);
    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 123);
    HIDE_OUTPUT([&] {
        command(fmt::format("rerun {} first 1 last 1 every 1 post no dump x y z", dump_file));
    });
    lmp->output->thermo->evaluate_keyword("pe", &pe_rerun);
    ASSERT_DOUBLE_EQ(pe_1, pe_rerun);
    HIDE_OUTPUT([&] {
        command(fmt::format("rerun {} first 2 last 2 every 1 post yes dump x y z", dump_file));
    });
    lmp->output->thermo->evaluate_keyword("pe", &pe_rerun);
    ASSERT_DOUBLE_EQ(pe_2, pe_rerun);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, rerun_bin)
{
    auto dump_file = binary_dump_filename("rerun");
    HIDE_OUTPUT([&] {
        command("fix 1 all nve");
    });
    generate_dump(dump_file, "", 1);
    double pe_1, pe_2, pe_rerun;
    lmp->output->thermo->evaluate_keyword("pe", &pe_1);
    ASSERT_FILE_EXISTS(dump_file);
    continue_dump(1);
    lmp->output->thermo->evaluate_keyword("pe", &pe_2);
    ASSERT_FILE_EXISTS(dump_file);
    HIDE_OUTPUT([&] {
        command(fmt::format("rerun {} first 1 last 1 every 1 post no dump x y z", dump_file));
    });
    lmp->output->thermo->evaluate_keyword("pe", &pe_rerun);
    ASSERT_NEAR(pe_1, pe_rerun,1.0e-14);
    HIDE_OUTPUT([&] {
        command(fmt::format("rerun {} first 2 last 2 every 1 post yes dump x y z", dump_file));
    });
    lmp->output->thermo->evaluate_keyword("pe", &pe_rerun);
    ASSERT_NEAR(pe_2, pe_rerun,1.0e-14);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, multi_file_run1)
{
    auto dump_file = dump_filename("run1_*");
    generate_dump(dump_file, "", 1);

    auto run1_0 = dump_filename("run1_0");
    auto run1_1 = dump_filename("run1_1");
    ASSERT_FILE_EXISTS(run1_0);
    ASSERT_FILE_EXISTS(run1_1);
    ASSERT_EQ(count_lines(run1_0), 41);
    ASSERT_EQ(count_lines(run1_1), 41);
    delete_file(run1_0);
    delete_file(run1_1);
}

TEST_F(DumpAtomTest, per_processor_file_run1)
{
    auto dump_file = dump_filename("run1_p%");
    generate_dump(dump_file, "", 1);

    auto run1_p0 = dump_filename("run1_p0");
    ASSERT_FILE_EXISTS(run1_p0);
    ASSERT_EQ(count_lines(run1_p0), 82);
    delete_file(run1_p0);
}

TEST_F(DumpAtomTest, per_processor_multi_file_run1)
{
    auto dump_file = dump_filename("run1_p%_*");
    generate_dump(dump_file, "", 1);

    auto run1_p0_0 = dump_filename("run1_p0_0");
    auto run1_p0_1 = dump_filename("run1_p0_1");
    ASSERT_FILE_EXISTS(run1_p0_0);
    ASSERT_FILE_EXISTS(run1_p0_1);
    ASSERT_EQ(count_lines(run1_p0_0), 41);
    ASSERT_EQ(count_lines(run1_p0_1), 41);
    delete_file(run1_p0_0);
    delete_file(run1_p0_1);
}

TEST_F(DumpAtomTest, dump_modify_scale_invalid)
{
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 1 dump.txt");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*Expected boolean parameter instead of 'xxx'.*",
                 command("dump_modify id scale xxx"););
}

TEST_F(DumpAtomTest, dump_modify_image_invalid)
{
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 1 dump.txt");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*Expected boolean parameter instead of 'xxx'.*",
                 command("dump_modify id image xxx"););
}

TEST_F(DumpAtomTest, dump_modify_invalid)
{
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 1 dump.txt");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*Illegal dump_modify command.*", command("dump_modify id true"););
}

TEST_F(DumpAtomTest, write_dump)
{
    auto reference = dump_filename("run0_ref");
    auto dump_file = fmt::format("write_{}", dump_filename("run*"));

    BEGIN_HIDE_OUTPUT();
    command(fmt::format("dump id all atom 1 {}", reference));
    command("dump_modify id scale no units yes");
    command("run 0");
    command(fmt::format("write_dump all atom {} modify scale no units yes", dump_file));
    END_HIDE_OUTPUT();

    dump_file = fmt::format("write_{}", dump_filename("run0"));
    ASSERT_FILE_EXISTS(reference);
    ASSERT_FILE_EXISTS(dump_file);

    ASSERT_FILE_EQUAL(reference, dump_file);
    delete_file(reference);
    delete_file(dump_file);
}

TEST_F(DumpAtomTest, binary_write_dump)
{
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto reference = binary_dump_filename("write_run0_ref");
    auto dump_file = fmt::format("write_{}", binary_dump_filename("write_dump_atom_run*_p%"));

    BEGIN_HIDE_OUTPUT();
    command(fmt::format("dump id all atom 1 {}", reference));
    command("dump_modify id scale no units yes");
    command("run 0");
    command(fmt::format("write_dump all atom {} modify scale no units yes", dump_file));
    END_HIDE_OUTPUT();

    dump_file = fmt::format("write_{}", binary_dump_filename("write_dump_atom_run0_p0"));
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
