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
#include "atom.h"
#include "domain.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "update.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>
#include <string>

using namespace LAMMPS_NS;

using testing::MatchesRegex;
using testing::StrEq;

using utils::sfgets;
using utils::sfread;
using utils::split_words;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

class FileOperationsTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "FileOperationsTest";
        LAMMPSTest::SetUp();
        ASSERT_NE(lmp, nullptr);

        FILE *fp = fopen("safe_file_read_test.txt", "wb");
        ASSERT_NE(fp, nullptr);
        fputs("one line\n", fp);
        fputs("two_lines\n", fp);
        fputs("\n", fp);
        fputs("no newline", fp);
        fclose(fp);
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        remove("safe_file_read_test.txt");
    }

    bool file_exists(const char *file)
    {
        FILE *fp = fopen(file, "r");
        fclose(fp);
        return fp ? true : false;
    }
};

#define MAX_BUF_SIZE 128
TEST_F(FileOperationsTest, safe_fgets)
{
    char buf[MAX_BUF_SIZE];

    FILE *fp = fopen("safe_file_read_test.txt", "r");
    ASSERT_NE(fp, nullptr);

    memset(buf, 0, MAX_BUF_SIZE);
    utils::sfgets(FLERR, buf, MAX_BUF_SIZE, fp, "safe_file_read_test.txt", lmp->error);
    ASSERT_THAT(buf, StrEq("one line\n"));

    memset(buf, 0, MAX_BUF_SIZE);
    utils::sfgets(FLERR, buf, MAX_BUF_SIZE, fp, "safe_file_read_test.txt", lmp->error);
    ASSERT_THAT(buf, StrEq("two_lines\n"));

    memset(buf, 0, MAX_BUF_SIZE);
    utils::sfgets(FLERR, buf, MAX_BUF_SIZE, fp, "safe_file_read_test.txt", lmp->error);
    ASSERT_THAT(buf, StrEq("\n"));

    memset(buf, 0, MAX_BUF_SIZE);
    utils::sfgets(FLERR, buf, 4, fp, "safe_file_read_test.txt", lmp->error);
    ASSERT_THAT(buf, StrEq("no "));

    memset(buf, 0, MAX_BUF_SIZE);
    utils::sfgets(FLERR, buf, MAX_BUF_SIZE, fp, "safe_file_read_test.txt", lmp->error);
    ASSERT_THAT(buf, StrEq("newline"));

    memset(buf, 0, MAX_BUF_SIZE);
    TEST_FAILURE(
        ".*ERROR on proc 0: Unexpected end of file while "
        "reading file 'safe_file_read_test.txt'.*",
        utils::sfgets(FLERR, buf, MAX_BUF_SIZE, fp, "safe_file_read_test.txt", lmp->error););
    fclose(fp);
}

#define MAX_BUF_SIZE 128
TEST_F(FileOperationsTest, safe_fread)
{
    char buf[MAX_BUF_SIZE];

    FILE *fp = fopen("safe_file_read_test.txt", "r");
    ASSERT_NE(fp, nullptr);

    memset(buf, 0, MAX_BUF_SIZE);
    utils::sfread(FLERR, buf, 1, 9, fp, "safe_file_read_test.txt", lmp->error);
    ASSERT_THAT(buf, StrEq("one line\n"));

    memset(buf, 0, MAX_BUF_SIZE);
    utils::sfread(FLERR, buf, 1, 10, fp, "safe_file_read_test.txt", nullptr);
    ASSERT_THAT(buf, StrEq("two_lines\n"));

    TEST_FAILURE(".*ERROR on proc 0: Unexpected end of file while "
                 "reading file 'safe_file_read_test.txt'.*",
                 utils::sfread(FLERR, buf, 1, 100, fp, "safe_file_read_test.txt", lmp->error););

    // short read but no error triggered due to passing a NULL pointer
    memset(buf, 0, MAX_BUF_SIZE);
    clearerr(fp);
    rewind(fp);
    utils::sfread(FLERR, buf, 1, 100, fp, "safe_file_read_test.txt", nullptr);
    ASSERT_THAT(buf, StrEq("one line\ntwo_lines\n\nno newline"));
    fclose(fp);
}

TEST_F(FileOperationsTest, logmesg)
{
    char buf[8];
    BEGIN_HIDE_OUTPUT();
    command("echo none");
    END_HIDE_OUTPUT();
    BEGIN_CAPTURE_OUTPUT();
    utils::logmesg(lmp, "one\n");
    command("log test_logmesg.log");
    utils::logmesg(lmp, "two\n");
    command("log none");
    std::string out = END_CAPTURE_OUTPUT();
    memset(buf, 0, 8);
    FILE *fp = fopen("test_logmesg.log", "r");
    fread(buf, 1, 8, fp);
    fclose(fp);
    ASSERT_THAT(out, StrEq("one\ntwo\n"));
    ASSERT_THAT(buf, StrEq("two\n"));
    remove("test_logmesg.log");
}

TEST_F(FileOperationsTest, restart)
{
    BEGIN_HIDE_OUTPUT();
    command("echo none");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Write_restart command before simulation box is defined.*",
                 command("write_restart test.restart"););

    BEGIN_HIDE_OUTPUT();
    command("region box block -2 2 -2 2 -2 2");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("mass 1 1.0");
    command("reset_timestep 333");
    command("comm_modify cutoff 0.2");
    command("write_restart noinit.restart noinit");
    command("run 0 post no");
    command("write_restart test.restart");
    command("write_restart step*.restart");
    command("write_restart multi-%.restart");
    command("write_restart multi2-%.restart fileper 2");
    command("write_restart multi3-%.restart nfile 1");
    if (info->has_package("MPIIO")) command("write_restart test.restart.mpiio");
    END_HIDE_OUTPUT();

    ASSERT_TRUE(file_exists("noinit.restart"));
    ASSERT_TRUE(file_exists("test.restart"));
    ASSERT_TRUE(file_exists("step333.restart"));
    ASSERT_TRUE(file_exists("multi-base.restart"));
    ASSERT_TRUE(file_exists("multi-0.restart"));
    ASSERT_TRUE(file_exists("multi2-base.restart"));
    ASSERT_TRUE(file_exists("multi2-0.restart"));
    ASSERT_TRUE(file_exists("multi3-base.restart"));
    ASSERT_TRUE(file_exists("multi3-0.restart"));
    if (info->has_package("MPIIO")) ASSERT_TRUE(file_exists("test.restart.mpiio"));

    if (!info->has_package("MPIIO")) {
        TEST_FAILURE(".*ERROR: Illegal write_restart command.*",
                     command("write_restart test.restart.mpiio"););
    } else {
        TEST_FAILURE(".*ERROR: Restart file MPI-IO output not allowed with % in filename.*",
                     command("write_restart test.restart-%.mpiio"););
    }

    TEST_FAILURE(".*ERROR: Illegal write_restart command.*", command("write_restart"););
    TEST_FAILURE(".*ERROR: Illegal write_restart command.*",
                 command("write_restart test.restart xxxx"););
    TEST_FAILURE(".*ERROR on proc 0: Cannot open restart file some_crazy_dir/test.restart:"
                 " No such file or directory.*",
                 command("write_restart some_crazy_dir/test.restart"););
    BEGIN_HIDE_OUTPUT();
    command("clear");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->update->ntimestep, 0);
    ASSERT_EQ(lmp->domain->triclinic, 0);

    TEST_FAILURE(
        ".*ERROR on proc 0: Cannot open restart file noexist.restart: No such file or directory.*",
        command("read_restart noexist.restart"););

    BEGIN_HIDE_OUTPUT();
    command("read_restart step333.restart");
    command("change_box all triclinic");
    command("write_restart triclinic.restart");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->atom->natoms, 1);
    ASSERT_EQ(lmp->update->ntimestep, 333);
    ASSERT_EQ(lmp->domain->triclinic, 1);
    BEGIN_HIDE_OUTPUT();
    command("clear");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->update->ntimestep, 0);
    ASSERT_EQ(lmp->domain->triclinic, 0);
    BEGIN_HIDE_OUTPUT();
    command("read_restart triclinic.restart");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->atom->natoms, 1);
    ASSERT_EQ(lmp->update->ntimestep, 333);
    ASSERT_EQ(lmp->domain->triclinic, 1);

    // clean up
    unlink("noinit.restart");
    unlink("test.restart");
    unlink("step333.restart");
    unlink("multi-base.restart");
    unlink("multi-0.restart");
    unlink("multi2-base.restart");
    unlink("multi2-0.restart");
    unlink("multi3-base.restart");
    unlink("multi3-0.restart");
    unlink("triclinic.restart");
    if (info->has_package("MPIIO")) unlink("test.restart.mpiio");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (Info::get_mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
