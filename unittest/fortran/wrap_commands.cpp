// unit tests for issuing command to a LAMMPS instance through the Fortran wrapper

#include "lammps.h"
#include <cstdio> // for stdin, stdout
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_file();
void f_lammps_command();
void f_lammps_commands_list();
void f_lammps_commands_string();
double f_lammps_get_natoms();
}

class LAMMPS_commands : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_commands(){};
    ~LAMMPS_commands() override{};

    void SetUp() override
    {
        ::testing::internal::CaptureStdout();
        lmp                = (LAMMPS_NS::LAMMPS *)f_lammps_with_args();
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 8).c_str(), "LAMMPS (");
    }
    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        f_lammps_close();
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};

TEST_F(LAMMPS_commands, from_file)
{
    EXPECT_EQ(f_lammps_get_natoms(), 0);
    f_lammps_file();
    EXPECT_EQ(f_lammps_get_natoms(), 2);
};

TEST_F(LAMMPS_commands, from_line)
{
    EXPECT_EQ(f_lammps_get_natoms(), 0);
    f_lammps_command();
    EXPECT_EQ(f_lammps_get_natoms(), 1);
};

TEST_F(LAMMPS_commands, from_list)
{
    EXPECT_EQ(f_lammps_get_natoms(), 0);
    f_lammps_commands_list();
    EXPECT_EQ(f_lammps_get_natoms(), 2);
};

TEST_F(LAMMPS_commands, from_string)
{
    EXPECT_EQ(f_lammps_get_natoms(), 0);
    f_lammps_commands_string();
    EXPECT_EQ(f_lammps_get_natoms(), 2);
};
