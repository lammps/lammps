// unit tests for the LAMMPS base class

#include "lammps.h"
#include <mpi.h>
#include <cstdio>  // for stdin, stdout
#include <string>

#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
    void *f_lammps_open_no_args();
    void *f_lammps_open_with_args();
    void f_lammps_close();
    int f_lammps_get_comm();
}

TEST(fortran_open, no_args) {
    ::testing::internal::CaptureStdout();
    void *handle = f_lammps_open_no_args();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0,6).c_str(),"LAMMPS");
    int mpi_init=0;
    MPI_Initialized(&mpi_init);
    EXPECT_GT(mpi_init,0);
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;
    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_NE(lmp->citeme, nullptr);
    ::testing::internal::CaptureStdout();
    f_lammps_close();
    output = testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0,16).c_str(), "Total wall time:");
}

TEST(fortran_open, with_args) {
    ::testing::internal::CaptureStdout();
    void *handle = f_lammps_open_with_args();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0,6).c_str(),"LAMMPS");
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);

    int f_comm = f_lammps_get_comm();
    MPI_Comm mycomm = MPI_Comm_f2c(f_comm);
    EXPECT_NE(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->world, mycomm);
    
    ::testing::internal::CaptureStdout();
    f_lammps_close();
    output = testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0,16).c_str(), "Total wall time:");
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
