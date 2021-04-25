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
#include "error.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "utils.h"
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
    char buf[64];
    BEGIN_HIDE_OUTPUT();
    command("echo none");
    END_HIDE_OUTPUT();
    BEGIN_CAPTURE_OUTPUT();
    utils::logmesg(lmp, "one\n");
    command("log test_logmesg.log");
    utils::logmesg(lmp, "two\n");
    utils::logmesg(lmp, "three={}\n", 3);
    utils::logmesg(lmp, "four {}\n");
    utils::logmesg(lmp, "five\n", 5);
    command("log none");
    std::string out = END_CAPTURE_OUTPUT();
    memset(buf, 0, 64);
    FILE *fp = fopen("test_logmesg.log", "r");
    fread(buf, 1, 64, fp);
    fclose(fp);
    ASSERT_THAT(out, StrEq("one\ntwo\nthree=3\nargument not found\nfive\n"));
    ASSERT_THAT(buf, StrEq("two\nthree=3\nargument not found\nfive\n"));
    remove("test_logmesg.log");
}

TEST_F(FileOperationsTest, error_message_warn)
{
    char buf[64];
    BEGIN_HIDE_OUTPUT();
    command("echo none");
    command("log test_error_warn.log");
    END_HIDE_OUTPUT();
    BEGIN_CAPTURE_OUTPUT();
    lmp->error->message("testme.cpp", 10, "message me");
    lmp->error->warning("testme.cpp", 100, "warn me");
    command("log none");
    std::string out = END_CAPTURE_OUTPUT();
    memset(buf, 0, 64);
    FILE *fp = fopen("test_error_warn.log", "r");
    fread(buf, 1, 64, fp);
    fclose(fp);
    auto msg = StrEq("message me (testme.cpp:10)\n"
                     "WARNING: warn me (testme.cpp:100)\n");
    ASSERT_THAT(out, msg);
    ASSERT_THAT(buf, msg);
    remove("test_error_warn.log");
}

TEST_F(FileOperationsTest, error_all_one)
{
    char buf[64];
    BEGIN_HIDE_OUTPUT();
    command("echo none");
    command("log none");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: exit \\(testme.cpp:10\\).*",
                 lmp->error->all("testme.cpp", 10, "exit"););
    TEST_FAILURE(".*ERROR: exit too \\(testme.cpp:10\\).*",
                 lmp->error->all("testme.cpp", 10, "exit {}", "too"););
    TEST_FAILURE(".*ERROR: argument not found \\(testme.cpp:10\\).*",
                 lmp->error->all("testme.cpp", 10, "exit {} {}", "too"););
    TEST_FAILURE(".*ERROR on proc 0: exit \\(testme.cpp:10\\).*",
                 lmp->error->one("testme.cpp", 10, "exit"););
    TEST_FAILURE(".*ERROR on proc 0: exit too \\(testme.cpp:10\\).*",
                 lmp->error->one("testme.cpp", 10, "exit {}", "too"););
    TEST_FAILURE(".*ERROR on proc 0: argument not found \\(testme.cpp:10\\).*",
                 lmp->error->one("testme.cpp", 10, "exit {} {}", "too"););
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
