// unit tests for issuing command to a LAMMPS instance through the library interface

#include "lammps.h"
#include "library.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

using ::testing::HasSubstr;
using ::testing::StartsWith;

const char *demo_input[] = {"region       box block 0 $x 0 2 0 2", "create_box 1 box",
                            "create_atoms 1 single 1.0 1.0 ${zpos}"};
const char *cont_input[] = {"create_atoms 1 single &", "0.2 0.1 0.1"};

class LibraryCommands : public ::testing::Test {
protected:
    void *lmp;
    LibraryCommands(){};
    ~LibraryCommands() override{};

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "screen", "-nocite",
                              "-var",        "x",    "2",    "-var",  "zpos",   "1.5"};

        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp                = lammps_open_no_mpi(argc, argv, NULL);
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        EXPECT_THAT(output, StartsWith("LAMMPS ("));
    }
    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        lammps_close(lmp);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, HasSubstr("Total wall time:"));
        if (verbose) std::cout << output;
        lmp = nullptr;
    }
};

TEST_F(LibraryCommands, from_file)
{
    FILE *fp;
    const char demo_file[] = "in.test";
    const char cont_file[] = "in.cont";

    fp = fopen(demo_file, "w");
    for (unsigned int i = 0; i < sizeof(demo_input) / sizeof(char *); ++i) {
        fputs(demo_input[i], fp);
        fputc('\n', fp);
    }
    fclose(fp);
    fp = fopen(cont_file, "w");
    for (unsigned int i = 0; i < sizeof(cont_input) / sizeof(char *); ++i) {
        fputs(cont_input[i], fp);
        fputc('\n', fp);
    }
    fclose(fp);

    EXPECT_EQ(lammps_get_natoms(lmp), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, demo_file);
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 1);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, cont_file);
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 2);

    unlink(demo_file);
    unlink(cont_file);
};

TEST_F(LibraryCommands, from_line)
{
    EXPECT_EQ(lammps_get_natoms(lmp), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    for (unsigned int i = 0; i < sizeof(demo_input) / sizeof(char *); ++i) {
        lammps_command(lmp, demo_input[i]);
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 1);
};

TEST_F(LibraryCommands, from_list)
{
    EXPECT_EQ(lammps_get_natoms(lmp), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_commands_list(lmp, sizeof(demo_input) / sizeof(char *), demo_input);
    lammps_commands_list(lmp, sizeof(cont_input) / sizeof(char *), cont_input);
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 2);
};

TEST_F(LibraryCommands, from_string)
{
    std::string cmds("");

    for (unsigned int i = 0; i < sizeof(demo_input) / sizeof(char *); ++i) {
        cmds += demo_input[i];
        cmds += "\n";
    }
    for (unsigned int i = 0; i < sizeof(cont_input) / sizeof(char *); ++i) {
        cmds += cont_input[i];
        cmds += "\n";
    }
    EXPECT_EQ(lammps_get_natoms(lmp), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_commands_string(lmp, cmds.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 2);

    // repeat test with DOS/Windows style CR-LF line endings
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "clear");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    cmds.clear();
    for (unsigned int i = 0; i < sizeof(demo_input) / sizeof(char *); ++i) {
        cmds += demo_input[i];
        cmds += "\r\n";
    }
    for (unsigned int i = 0; i < sizeof(cont_input) / sizeof(char *); ++i) {
        cmds += cont_input[i];
        cmds += "\r\n";
    }
    EXPECT_EQ(lammps_get_natoms(lmp), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_commands_string(lmp, cmds.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 2);
};
