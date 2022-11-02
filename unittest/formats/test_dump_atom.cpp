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

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include "../testing/utils.h"
#include "fmt/format.h"
#include "output.h"
#include "thermo.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <string>
#include <vector>

using ::testing::Eq;

char *BINARY2TXT_EXECUTABLE = nullptr;
bool verbose                = false;

namespace LAMMPS_NS {

class DumpAtomTest : public MeltTest {
    std::string dump_style = "atom";

public:
    void enable_triclinic()
    {
        BEGIN_HIDE_OUTPUT();
        command("change_box all triclinic");
        END_HIDE_OUTPUT();
    }

    std::string dump_filename(const std::string &ident)
    {
        return fmt::format("dump_{}_{}.melt", dump_style, ident);
    }

    std::string text_dump_filename(const std::string &ident)
    {
        return fmt::format("dump_{}_text_{}.melt", dump_style, ident);
    }

    std::string binary_dump_filename(const std::string &ident)
    {
        return fmt::format("dump_{}_binary_{}.melt.bin", dump_style, ident);
    }

    void generate_dump(const std::string &dump_file, const std::string &dump_modify_options,
                       int ntimesteps)
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

    void close_dump()
    {
        BEGIN_HIDE_OUTPUT();
        command("undump id");
        END_HIDE_OUTPUT();
    }

    void generate_text_and_binary_dump(const std::string &text_file, const std::string &binary_file,
                                       const std::string &dump_modify_options, int ntimesteps)
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

    std::string convert_binary_to_text(const std::string &binary_file)
    {
        BEGIN_HIDE_OUTPUT();
        std::string cmdline = fmt::format("\"{}\" {}", BINARY2TXT_EXECUTABLE, binary_file);
        system(cmdline.c_str());
        END_HIDE_OUTPUT();
        return fmt::format("{}.txt", binary_file);
    }

    std::vector<std::string> extract_items(const std::string &file, const std::string &item)
    {
        std::string match = fmt::format("^ITEM: {}$", item);
        std::vector<std::string> values;

        std::ifstream dump(file);
        for (std::string buffer; std::getline(dump, buffer); /* */) {
            buffer = utils::trim(buffer);
            if (utils::strmatch(buffer, match)) {
                std::getline(dump, buffer);
                values.push_back(buffer);
            }
        }
        return values;
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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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
    close_dump();
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
    close_dump();
    lmp->output->thermo->evaluate_keyword("pe", &pe_2);
    ASSERT_FILE_EXISTS(dump_file);
    HIDE_OUTPUT([&] {
        command(fmt::format("rerun {} first 1 last 1 every 1 post no dump x y z", dump_file));
    });
    lmp->output->thermo->evaluate_keyword("pe", &pe_rerun);
    ASSERT_NEAR(pe_1, pe_rerun, 1.0e-14);
    HIDE_OUTPUT([&] {
        command(fmt::format("rerun {} first 2 last 2 every 1 post yes dump x y z", dump_file));
    });
    lmp->output->thermo->evaluate_keyword("pe", &pe_rerun);
    ASSERT_NEAR(pe_2, pe_rerun, 1.0e-14);
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

    TEST_FAILURE(".*Unknown dump_modify keyword: true.*", command("dump_modify id true"););
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
    if (!BINARY2TXT_EXECUTABLE) GTEST_SKIP();

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

TEST_F(DumpAtomTest, frequency)
{
    auto dump_file = dump_filename("frequency");
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 5 " + dump_file);
    command("run 15 post no");
    command("run 12 post no");
    END_HIDE_OUTPUT();

    // NOTE: must reset to current timestep (27) to avoid unexpected issues with following
    TEST_FAILURE(".*ERROR: Cannot reset timestep with active dump - must undump first.*",
                 command("reset_timestep 27"););

    BEGIN_HIDE_OUTPUT();
    command("run 3 post no");
    command("undump id");
    command("reset_timestep 5");
    command("dump id all atom 10 " + dump_file);
    command("dump_modify id append yes");
    command("run 20 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    std::vector<std::string> expected, values;
    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"0", "5", "10", "15", "20", "25", "30", "10", "20"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    BEGIN_HIDE_OUTPUT();
    command("reset_timestep 10");
    command("dump id all atom 10 " + dump_file);
    command("run 20 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"10", "20", "30"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    BEGIN_HIDE_OUTPUT();
    command("reset_timestep 0");
    command("dump id all atom 10 " + dump_file);
    command("minimize 0.0 0.0 15 30");
    command("run 20 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"0", "10", "15", "20", "30"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    delete_file(dump_file);
}

//-------------------------------------------------------------------------------------------------
// dump_modify
//-------------------------------------------------------------------------------------------------

TEST_F(DumpAtomTest, delay)
{
    auto dump_file = dump_filename("delay");
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 10 " + dump_file);
    command("dump_modify id delay 20");
    command("run 50 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    std::vector<std::string> expected, values;
    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"20", "30", "40", "50"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    delete_file(dump_file);
}

TEST_F(DumpAtomTest, colname)
{
    auto dump_file = dump_filename("colname");
    BEGIN_HIDE_OUTPUT();
    command("group one id 1");
    command("dump id one atom 10 " + dump_file);
    command("run 5 post no");
    command("dump_modify id colname id AtomID colname 3 x-scaled colname -4 z-scaled");
    command("run 10 post no");
    command("dump_modify id colname default");
    command("run 10 post no");
    command("dump_modify id colname id AtomID colname 3 x-scaled colname -4 z-scaled");
    command("dump_modify id scale no image yes");
    command("run 10 post no");
    command("dump_modify id colname id AtomID colname 3 X colname -4 Z colname ix img_x");
    command("run 10 post no");
    command("dump_modify id colname default");
    command("run 10 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    std::vector<std::string> expected, values;
    values   = extract_items(dump_file, "ATOMS id type xs ys zs");
    expected = {"1 1 0 0 0", "1 1 0 0 0"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    values   = extract_items(dump_file, "ATOMS AtomID type x-scaled ys z-scaled");
    expected = {"1 1 0 0 0"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    values   = extract_items(dump_file, "ATOMS id type x y z ix iy iz");
    expected = {"1 1 0 0 0 0 0 0", "1 1 0 0 0 0 0 0"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    values   = extract_items(dump_file, "ATOMS AtomID type X y Z img_x iy iz");
    expected = {"1 1 0 0 0 0 0 0"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    delete_file(dump_file);
}

TEST_F(DumpAtomTest, units_time)
{
    auto dump_file = dump_filename("units_time");
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 10 " + dump_file);
    command("dump_modify id units yes time yes");
    command("run 30 post no");
    command("timestep 0.01");
    command("run 30 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    std::vector<std::string> expected, values;
    values   = extract_items(dump_file, "TIME");
    expected = {"0", "0.05", "0.1", "0.15", "0.25", "0.35", "0.45"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    values   = extract_items(dump_file, "UNITS");
    expected = {"lj"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    delete_file(dump_file);
}

TEST_F(DumpAtomTest, every)
{
    auto dump_file = dump_filename("every");
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 10 " + dump_file);
    command("run 20 post no");
    command("dump_modify id every 5");
    command("run 15 post no");
    command("dump_modify id every 10");
    command("run 25 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    std::vector<std::string> expected, values;
    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"0", "10", "20", "25", "30", "35", "40", "50", "60"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    delete_file(dump_file);

    BEGIN_HIDE_OUTPUT();
    command("reset_timestep 0");
    command("dump id all atom 1 " + dump_file);
    command("variable next equal (step+1)*(step+1)");
    command("dump_modify id every v_next");
    command("run 50 post no");
    command("variable next equal logfreq(10,7,10)");
    command("dump_modify id every v_next");
    command("run 100 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"1", "4", "25", "60", "70", "100"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    delete_file(dump_file);
}

TEST_F(DumpAtomTest, every_time)
{
    auto dump_file = dump_filename("every_time");
    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 10 " + dump_file);
    command("dump_modify id every/time 0.1");
    command("run 40 post no");
    command("timestep 0.01");
    command("run 20 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    std::vector<std::string> expected, values;
    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"0", "20", "40", "50", "60"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    delete_file(dump_file);
}

TEST_F(DumpAtomTest, header)
{
    auto dump_file = dump_filename("header");
    BEGIN_HIDE_OUTPUT();
    command("reset_timestep 5");
    command("dump id all atom 10 " + dump_file);
    command("dump_modify id first no header yes");
    command("run 40 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    std::vector<std::string> expected, values;
    values   = extract_items(dump_file, "TIMESTEP");
    expected = {"10", "20", "30", "40"};
    ASSERT_EQ(values.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        ASSERT_THAT(values[i], Eq(expected[i]));

    BEGIN_HIDE_OUTPUT();
    command("dump id all atom 10 " + dump_file);
    command("dump_modify id header no");
    command("run 40 post no");
    command("undump id");
    END_HIDE_OUTPUT();

    values = extract_items(dump_file, "TIMESTEP");
    ASSERT_EQ(values.size(), 0);
    delete_file(dump_file);
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    BINARY2TXT_EXECUTABLE = getenv("BINARY2TXT_EXECUTABLE");

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
