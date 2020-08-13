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

using ::testing::MatchesRegex;

#define TEST_FAILURE(errmsg, ...)                                 \
    if (Info::has_exceptions()) {                                 \
        ::testing::internal::CaptureStdout();                     \
        ASSERT_ANY_THROW({__VA_ARGS__});                          \
        auto mesg = ::testing::internal::GetCapturedStdout();     \
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (Info::get_mpi_vendor() != "Open MPI") {               \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

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

class DumpCustomTest : public MeltTest {
};

TEST_F(DumpCustomTest, run0)
{
    auto dump_file = "dump_custom_run0.melt";

    if (!verbose) ::testing::internal::CaptureStdout();
    command(fmt::format("dump id all custom 1 {} id type x y vx fx", dump_file));
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();


    ASSERT_FILE_EXISTS(dump_file);
    auto lines = read_lines(dump_file);
    ASSERT_EQ(lines.size(), 41);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 2);
    ASSERT_STREQ(lines[8].c_str(), "ITEM: ATOMS id type x y vx fx");
    ASSERT_EQ(utils::split_words(lines[9]).size(), 6);
    delete_file(dump_file);
}

TEST_F(DumpCustomTest, triclinic_run0)
{
    auto dump_file = "dump_custom_tri_run0.melt";
    if (!verbose) ::testing::internal::CaptureStdout();

    command("change_box all triclinic");
    command(fmt::format("dump id all custom 1 {} id type x y vx fx", dump_file));
    command("run 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FILE_EXISTS(dump_file);

    auto lines = read_lines(dump_file);
    ASSERT_STREQ(lines[4].c_str(), "ITEM: BOX BOUNDS xy xz yz pp pp pp");
    ASSERT_EQ(utils::split_words(lines[5]).size(), 3);

    ASSERT_EQ(lines.size(), 41);
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
