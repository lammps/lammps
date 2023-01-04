// unit tests for the LAMMPS base class

#include "lammps.h"
#include <cstdio> // for stdin, stdout
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
void *f_lammps_open_no_args();
void *f_lammps_open_with_args();
void *f_lammps_no_mpi_no_args();
void *f_lammps_no_mpi_with_args();
void f_lammps_close();
int f_lammps_get_comm();
}

// C wrapper to split MPI communicator w/o requiring a Fortran MPI lib
extern "C" int create_mpi_comm_split(int color, int key)
{
    MPI_Comm c_newcomm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &c_newcomm);
    return MPI_Comm_c2f(c_newcomm);
}

TEST(open_no_mpi, no_args)
{
    ::testing::internal::CaptureStdout();
    int mpi_init = 0;
    MPI_Initialized(&mpi_init);
    EXPECT_EQ(mpi_init, 0);
    void *handle       = f_lammps_no_mpi_no_args();
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 6).c_str(), "LAMMPS");
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;
    MPI_Initialized(&mpi_init);
    EXPECT_NE(mpi_init, 0);
    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_NE(lmp->citeme, nullptr);
    ::testing::internal::CaptureStdout();
    f_lammps_close();
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
}

TEST(open_no_mpi, with_args)
{
    ::testing::internal::CaptureStdout();
    void *handle       = f_lammps_no_mpi_with_args();
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 6).c_str(), "LAMMPS");
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);
    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);

    ::testing::internal::CaptureStdout();
    f_lammps_close();
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
}

TEST(fortran_open, no_args)
{
    ::testing::internal::CaptureStdout();
    void *handle       = f_lammps_open_no_args();
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 6).c_str(), "LAMMPS");
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    int f_comm      = f_lammps_get_comm();
    MPI_Comm mycomm = MPI_Comm_f2c(f_comm);
    EXPECT_EQ(lmp->world, mycomm);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_NE(lmp->citeme, nullptr);
    ::testing::internal::CaptureStdout();
    f_lammps_close();
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
}

TEST(fortran_open, with_args)
{
    ::testing::internal::CaptureStdout();
    void *handle       = f_lammps_open_with_args();
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 6).c_str(), "LAMMPS");
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    int f_comm      = f_lammps_get_comm();
    MPI_Comm mycomm = MPI_Comm_f2c(f_comm);
    EXPECT_EQ(lmp->world, mycomm);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);

    ::testing::internal::CaptureStdout();
    f_lammps_close();
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
}
