// unit tests for checking LAMMPS configuration settings  through the library interface

#include "lammps.h"
#include "library.h"
#include "timer.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/test_mpi_main.h"

using ::testing::ExitedWithCode;
using ::testing::HasSubstr;
using ::testing::StartsWith;
using ::testing::StrEq;

TEST(MPI, global_box)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    EXPECT_EQ(nprocs, 4);
    EXPECT_GT(me, -1);
    EXPECT_LT(me, 5);

    double boxlo[3];
    double boxhi[3];
    double xy  = 0.0;
    double yz  = 0.0;
    double xz  = 0.0;
    int pflags[3];
    int boxflag;

    ::testing::internal::CaptureStdout();
    const char *args[] = {"LAMMPS_test", "-log",      "none",
                          "-echo",       "screen",    "-nocite"};
    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    void * lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
    lammps_command(lmp, "units           lj");
    lammps_command(lmp, "atom_style      atomic");
    lammps_command(lmp, "region          box block 0 2 0 2 0 2");
    lammps_command(lmp, "create_box      1 box");

    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(boxlo[0], 0.0);
    EXPECT_EQ(boxlo[1], 0.0);
    EXPECT_EQ(boxlo[2], 0.0);

    EXPECT_EQ(boxhi[0], 2.0);
    EXPECT_EQ(boxhi[1], 2.0);
    EXPECT_EQ(boxhi[2], 2.0);

    ::testing::internal::CaptureStdout();
    lammps_close(lmp);
    ::testing::internal::GetCapturedStdout();
};

TEST(MPI, sub_box)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    EXPECT_EQ(nprocs, 4);
    EXPECT_GT(me, -1);
    EXPECT_LT(me, 5);

    double boxlo[3];
    double boxhi[3];
    double xy  = 0.0;
    double yz  = 0.0;
    double xz  = 0.0;
    int pflags[3];
    int boxflag;

    ::testing::internal::CaptureStdout();
    const char *args[] = {"LAMMPS_test", "-log",      "none",
                          "-echo",       "screen",    "-nocite"};
    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    void * lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
    lammps_command(lmp, "units           lj");
    lammps_command(lmp, "atom_style      atomic");
    lammps_command(lmp, "region          box block 0 2 0 2 0 2");
    lammps_command(lmp, "create_box      1 box");

    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(boxlo[0], 0.0);
    EXPECT_EQ(boxlo[1], 0.0);
    EXPECT_EQ(boxlo[2], 0.0);

    EXPECT_EQ(boxhi[0], 2.0);
    EXPECT_EQ(boxhi[1], 2.0);
    EXPECT_EQ(boxhi[2], 2.0);

    double * sublo = (double*)lammps_extract_global(lmp, "sublo");
    double * subhi = (double*)lammps_extract_global(lmp, "subhi");

    ASSERT_NE(sublo, nullptr);
    ASSERT_NE(subhi, nullptr);

    EXPECT_GE(sublo[0], boxlo[0]);
    EXPECT_GE(sublo[1], boxlo[1]);
    EXPECT_GE(sublo[2], boxlo[2]);
    EXPECT_LE(subhi[0], boxhi[0]);
    EXPECT_LE(subhi[1], boxhi[1]);
    EXPECT_LE(subhi[2], boxhi[2]);

    ::testing::internal::CaptureStdout();
    lammps_close(lmp);
    ::testing::internal::GetCapturedStdout();
};
