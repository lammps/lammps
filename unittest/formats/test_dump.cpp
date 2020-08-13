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

#include "info.h"
#include "input.h"
#include "lammps.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "fmt/format.h"
#include "potential_file_reader.h"

#include <cstring>
#include <mpi.h>
#include <iostream>
#include <fstream>

using namespace LAMMPS_NS;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

char * BINARY2TXT_BINARY = nullptr;

static void delete_file(const std::string &filename)
{
    remove(filename.c_str());
}

static size_t count_lines(const std::string &filename)
{
    std::ifstream infile(filename);
    std::string line;
    size_t nlines = 0;

    while (std::getline(infile, line))
        ++nlines;

    return nlines;
}

static bool equal_lines(const std::string &fileA, const std::string &fileB)
{
    std::ifstream afile(fileA);
    std::ifstream bfile(fileB);
    std::string lineA, lineB;

    while (std::getline(afile, lineA)) {
        if(!std::getline(bfile, lineB)) return false;
        if(lineA != lineB) return false;
    }

    return true;
}

static std::vector<std::string> read_lines(const std::string &filename) {
    std::vector<std::string> lines;
    std::ifstream infile(filename);
    std::string line;

    while (std::getline(infile, line))
        lines.push_back(line);

    return lines;
}

static bool file_exists(const std::string &filename) {
    struct stat result;
    return stat(filename.c_str(), &result) == 0;
}

#define ASSERT_FILE_EXISTS(NAME) ASSERT_TRUE(file_exists(NAME))
#define ASSERT_FILE_EQUAL(FILE_A, FILE_B) ASSERT_TRUE(equal_lines(FILE_A, FILE_B))

class LAMMPSTest : public ::testing::Test {
public:
    void command(const std::string &line) {
        lmp->input->one(line.c_str());
    }

protected:
    const char * testbinary = "LAMMPSTest";
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = { testbinary, "-log", "none", "-echo", "screen", "-nocite"};
        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        InitSystem();
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }


    virtual void InitSystem() {
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

class MeltTest : public LAMMPSTest {
protected:
    virtual void InitSystem() override {
        command("units           lj");
        command("atom_style      atomic");

        command("lattice         fcc 0.8442");
        command("region          box block 0 2 0 2 0 2");
        command("create_box      1 box");
        command("create_atoms    1 box");
        command("mass            1 1.0");

        command("velocity        all create 3.0 87287");

        command("pair_style      lj/cut 2.5");
        command("pair_coeff      1 1 1.0 1.0 2.5");

        command("neighbor        0.3 bin");
        command("neigh_modify    every 20 delay 0 check no");        
    }
};

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

TEST_F(DumpAtomTest, tricilinic_run0)
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
