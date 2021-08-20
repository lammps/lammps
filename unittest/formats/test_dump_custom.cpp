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

using ::testing::Eq;

char *BINARY2TXT_BINARY = nullptr;
bool verbose            = false;

class DumpCustomTest : public MeltTest {
    std::string dump_style = "custom";

public:
    void enable_triclinic()
    {
        BEGIN_HIDE_OUTPUT();
        command("change_box all triclinic");
        END_HIDE_OUTPUT();
    }

    void generate_dump(std::string dump_file, std::string fields, std::string dump_modify_options,
                       int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id all {} 1 {} {}", dump_style, dump_file, fields));

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
                                       std::string fields, std::string dump_modify_options,
                                       int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id0 all {} 1 {} {}", dump_style, text_file, fields));
        command(fmt::format("dump id1 all {} 1 {} {}", dump_style, binary_file, fields));

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

TEST_F(DumpCustomTest, run1)
{
    auto dump_file = "dump_custom_run1.melt";
    auto fields =
        "id type proc procp1 mass x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    generate_dump(dump_file, fields, "units yes", 1);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 84);
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[10], Eq(fmt::format("ITEM: ATOMS {}", fields)));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 26);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, thresh_run0)
{
    auto dump_file = "dump_custom_thresh_run0.melt";
    auto fields    = "id type x y z";

    generate_dump(dump_file, fields, "units yes thresh x < 1 thresh y < 1 thresh z < 1", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 15);
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[10], Eq(fmt::format("ITEM: ATOMS {}", fields)));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 5);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, compute_run0)
{
    BEGIN_HIDE_OUTPUT();
    command("compute comp all property/atom x y z");
    END_HIDE_OUTPUT();

    auto dump_file = "dump_custom_compute_run0.melt";
    auto fields    = "id type x y z c_comp[1] c_comp[2] c_comp[3]";

    generate_dump(dump_file, fields, "units yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 43);
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[10], Eq(fmt::format("ITEM: ATOMS {}", fields)));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 8);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, fix_run0)
{
    if (!info->has_style("fix", "numdiff")) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("fix numdiff all numdiff 1 0.0001");
    END_HIDE_OUTPUT();

    auto dump_file = "dump_custom_compute_run0.melt";
    auto fields    = "id x y z f_numdiff[1] f_numdiff[2] f_numdiff[3]";

    generate_dump(dump_file, fields, "units yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 43);
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[10], Eq(fmt::format("ITEM: ATOMS {}", fields)));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 7);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, custom_run0)
{
    BEGIN_HIDE_OUTPUT();
    command("fix prop all property/atom i_flag1 d_flag2");
    command("compute 1 all property/atom i_flag1 d_flag2");
    END_HIDE_OUTPUT();

    auto dump_file = "dump_custom_custom_run0.melt";
    auto fields    = "id x y z i_flag1 d_flag2";

    generate_dump(dump_file, fields, "units yes", 0);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 43);
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[10], Eq(fmt::format("ITEM: ATOMS {}", fields)));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 6);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, binary_run1)
{
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = "dump_custom_text_run1.melt";
    auto binary_file = "dump_custom_binary_run1.melt.bin";
    auto fields = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    generate_text_and_binary_dump(text_file, binary_file, fields, "units yes", 1);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpCustomTest, triclinic_run1)
{
    auto dump_file = "dump_custom_tri_run1.melt";
    auto fields    = "id type proc x y z ix iy iz xs ys zs xu yu zu xsu ysu zsu vx vy vz fx fy fz";

    enable_triclinic();

    generate_dump(dump_file, fields, "units yes", 1);

    ASSERT_FILE_EXISTS(dump_file);

    auto lines = read_lines(dump_file);
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS xy xz yz pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 3);

    ASSERT_EQ(lines.size(), 84);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, binary_triclinic_run1)
{
    if (!BINARY2TXT_BINARY) GTEST_SKIP();

    auto text_file   = "dump_custom_tri_text_run1.melt";
    auto binary_file = "dump_custom_tri_binary_run1.melt.bin";
    auto fields      = "id type proc x y z xs ys zs xsu ysu zsu vx vy vz fx fy fz";

    enable_triclinic();

    generate_text_and_binary_dump(text_file, binary_file, fields, "units yes", 1);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(binary_file);

    auto converted_file = convert_binary_to_text(binary_file);

    ASSERT_FILE_EXISTS(converted_file);
    ASSERT_FILE_EQUAL(text_file, converted_file);
    delete_file(text_file);
    delete_file(binary_file);
    delete_file(converted_file);
}

TEST_F(DumpCustomTest, with_variable_run1)
{
    BEGIN_HIDE_OUTPUT();
    command("compute         1 all property/atom proc");
    command("variable        p atom (c_1%10)+1");
    END_HIDE_OUTPUT();

    auto dump_file = "dump_custom_with_variable_run1.melt";
    auto fields    = "id type x y z v_p";

    generate_dump(dump_file, fields, "units yes", 1);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 84);
    ASSERT_THAT(lines[6], Eq("ITEM: BOX BOUNDS pp pp pp"));
    ASSERT_EQ(utils::split_words(lines[7]).size(), 2);
    ASSERT_THAT(lines[10], Eq("ITEM: ATOMS id type x y z v_p"));
    ASSERT_EQ(utils::split_words(lines[11]).size(), 6);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, run1plus1)
{
    auto dump_file = "dump_custom_run1plus1.melt";
    auto fields    = "id type x y z";

    generate_dump(dump_file, fields, "units yes", 1);

    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 84);
    continue_dump(1);
    ASSERT_FILE_EXISTS(dump_file);
    lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 125);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, run2)
{
    auto dump_file = "dump_custom_run2.melt";
    auto fields    = "id type x y z";
    generate_dump(dump_file, fields, "", 2);

    ASSERT_FILE_EXISTS(dump_file);
    ASSERT_EQ(count_lines(dump_file), 123);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, rerun)
{
    auto dump_file = "dump_rerun.melt";
    auto fields    = "id type xs ys zs";

    HIDE_OUTPUT([&] {
        command("fix 1 all nve");
    });
    generate_dump(dump_file, fields, "format float %20.15g", 1);
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
