// unit tests for issuing command to a LAMMPS instance through the Fortran wrapper

#include "lammps.h"
#include <mpi.h>
#include <cstdio>  // for stdin, stdout
#include <string>

#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
    void *f_lammps_with_args();
    void f_lammps_close();
    void f_lammps_file();
    double f_lammps_get_natoms();
}

class LAMMPS_commands : public ::testing::Test
{
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_commands() {};
    ~LAMMPS_commands() override {};

    void SetUp() override {
        ::testing::internal::CaptureStdout();
        lmp = (LAMMPS_NS::LAMMPS *)f_lammps_with_args();
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0,8).c_str(), "LAMMPS (");
    }
    void TearDown() override {
        ::testing::internal::CaptureStdout();
        f_lammps_close();
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0,16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};            

TEST_F(LAMMPS_commands, from_file) {
    EXPECT_EQ(f_lammps_get_natoms(),0);
    f_lammps_file();
    EXPECT_EQ(f_lammps_get_natoms(),2);
};

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
