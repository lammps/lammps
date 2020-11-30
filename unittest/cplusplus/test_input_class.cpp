// unit tests for issuing command to a LAMMPS instance through the Input class

#include "atom.h"
#include "input.h"
#include "lammps.h"
#include "memory.h"
#include <cstring>
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

const char *demo_input[] = {"region       box block 0 $x 0 2 0 2", "create_box 1 box",
                            "create_atoms 1 single 1.0 1.0 ${zpos}"};
const char *cont_input[] = {"create_atoms 1 single &", "0.2 0.1 0.1"};

namespace LAMMPS_NS {

class Input_commands : public ::testing::Test {
protected:
    LAMMPS *lmp;
    Input_commands()
    {
        const char *args[] = {"LAMMPS_test"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        int flag;
        MPI_Initialized(&flag);
        if (!flag) MPI_Init(&argc, &argv);
    }
    ~Input_commands() override {}

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "screen", "-nocite",
                              "-var",        "zpos", "1.5",  "-var",  "x",      "2"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp                = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 8).c_str(), "LAMMPS (");
    }
    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        delete lmp;
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};

TEST_F(Input_commands, from_file)
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

    EXPECT_EQ(lmp->atom->natoms, 0);
    lmp->input->file(demo_file);
    lmp->input->file(cont_file);
    EXPECT_EQ(lmp->atom->natoms, 2);

    unlink(demo_file);
    unlink(cont_file);
};

TEST_F(Input_commands, from_line)
{
    EXPECT_EQ(lmp->atom->natoms, 0);
    for (unsigned int i = 0; i < sizeof(demo_input) / sizeof(char *); ++i) {
        lmp->input->one(demo_input[i]);
    }
    EXPECT_EQ(lmp->atom->natoms, 1);
};

TEST_F(Input_commands, substitute)
{
    char *string, *scratch;
    int nstring = 100, nscratch = 100;

    lmp->memory->create(string, nstring, "test:string");
    lmp->memory->create(scratch, nscratch, "test:scratch");
    strcpy(string, demo_input[0]);
    lmp->input->substitute(string, scratch, nstring, nscratch, 0);
    EXPECT_STREQ(string, "region       box block 0 2 0 2 0 2");

    strcpy(string, demo_input[2]);
    lmp->input->substitute(string, scratch, nstring, nscratch, 0);
    EXPECT_STREQ(string, "create_atoms 1 single 1.0 1.0 1.5");
    lmp->memory->destroy(string);
    lmp->memory->destroy(scratch);
};
} // namespace LAMMPS_NS
