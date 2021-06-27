// unit tests for checking LAMMPS configuration settings  through the library interface

#define LAMMPS_LIB_MPI 1
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
    double xy = 0.0;
    double yz = 0.0;
    double xz = 0.0;
    int pflags[3];
    int boxflag;

    ::testing::internal::CaptureStdout();
    const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "screen", "-nocite"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);
    void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
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
    double xy = 0.0;
    double yz = 0.0;
    double xz = 0.0;
    int pflags[3];
    int boxflag;

    ::testing::internal::CaptureStdout();
    const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "screen", "-nocite"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);
    void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
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

    double *sublo = (double *)lammps_extract_global(lmp, "sublo");
    double *subhi = (double *)lammps_extract_global(lmp, "subhi");

    ASSERT_NE(sublo, nullptr);
    ASSERT_NE(subhi, nullptr);

    EXPECT_GE(sublo[0], boxlo[0]);
    EXPECT_GE(sublo[1], boxlo[1]);
    EXPECT_GE(sublo[2], boxlo[2]);
    EXPECT_LE(subhi[0], boxhi[0]);
    EXPECT_LE(subhi[1], boxhi[1]);
    EXPECT_LE(subhi[2], boxhi[2]);

    ::testing::internal::CaptureStdout();
    lammps_command(lmp, "change_box all triclinic");
    ::testing::internal::GetCapturedStdout();

    sublo = (double *)lammps_extract_global(lmp, "sublo_lambda");
    subhi = (double *)lammps_extract_global(lmp, "subhi_lambda");
    ASSERT_NE(sublo, nullptr);
    ASSERT_NE(subhi, nullptr);

    EXPECT_GE(sublo[0], 0.0);
    EXPECT_GE(sublo[1], 0.0);
    EXPECT_GE(sublo[2], 0.0);
    EXPECT_LE(subhi[0], 1.0);
    EXPECT_LE(subhi[1], 1.0);
    EXPECT_LE(subhi[2], 1.0);

    ::testing::internal::CaptureStdout();
    lammps_close(lmp);
    ::testing::internal::GetCapturedStdout();
};

TEST(MPI, split_comm)
{
    int nprocs, me, color, key;
    MPI_Comm newcomm;
    lammps_mpi_init();
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    color = me % 2;
    key   = me;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &newcomm);

    const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "screen", "-nocite"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);
    void *lmp          = lammps_open(argc, argv, newcomm, nullptr);
    lammps_command(lmp, "units           lj");
    lammps_command(lmp, "atom_style      atomic");
    lammps_command(lmp, "region          box block 0 2 0 2 0 2");
    lammps_command(lmp, "create_box      1 box");

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &me);
    EXPECT_EQ(nprocs, 2);
    EXPECT_GT(me, -1);
    EXPECT_LT(me, 2);
    EXPECT_EQ(lammps_extract_setting(lmp, "universe_size"), nprocs);
    EXPECT_EQ(lammps_extract_setting(lmp, "universe_rank"), me);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), nprocs);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_rank"), me);

    lammps_close(lmp);
};

TEST(MPI, multi_partition)
{
    FILE *fp;
    int nprocs, me;
    lammps_mpi_init();
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    const char *args[] = {"LAMMPS_test", "-log",   "none",    "-partition", "4x1",
                          "-echo",       "screen", "-nocite", "-in",        "none"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);
    void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);

    lammps_command(lmp, "units           lj");
    lammps_command(lmp, "atom_style      atomic");
    lammps_command(lmp, "region          box block 0 2 0 2 0 2");
    lammps_command(lmp, "create_box      1 box");
    lammps_command(lmp, "variable        partition universe 1 2 3 4");

    EXPECT_EQ(lammps_extract_setting(lmp, "universe_size"), nprocs);
    EXPECT_EQ(lammps_extract_setting(lmp, "universe_rank"), me);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_rank"), 0);

    char *part_id = (char *)lammps_extract_variable(lmp, "partition", nullptr);
    ASSERT_THAT(part_id, StrEq(std::to_string(me + 1)));

    lammps_close(lmp);
};

class MPITest : public ::testing::Test {
public:
    void command(const std::string &line) { lammps_command(lmp, line.c_str()); }

protected:
    const char *testbinary = "LAMMPSTest";
    void *lmp;

    void SetUp() override
    {
        const char *args[] = {testbinary, "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
        InitSystem();
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    virtual void InitSystem()
    {
        command("units           lj");
        command("atom_style      atomic");
        command("atom_modify     map yes");

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

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        lammps_close(lmp);
        lmp = nullptr;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(MPITest, size_rank)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    EXPECT_EQ(nprocs, lammps_extract_setting(lmp, "world_size"));
    EXPECT_EQ(me, lammps_extract_setting(lmp, "world_rank"));
}

#if !defined(LAMMPS_BIGBIG)

TEST_F(MPITest, gather)
{
    int64_t natoms = (int64_t)lammps_get_natoms(lmp);
    ASSERT_EQ(natoms, 32);
    int *p_nlocal = (int *)lammps_extract_global(lmp, "nlocal");
    int nlocal    = *p_nlocal;
    EXPECT_LT(nlocal, 32);
    EXPECT_EQ(nlocal, 8);

    // get the entire x on all procs
    double *x = new double[natoms * 3];
    lammps_gather(lmp, (char *)"x", 1, 3, x);

    int *tag         = (int *)lammps_extract_atom(lmp, "id");
    double **x_local = (double **)lammps_extract_atom(lmp, "x");

    // each proc checks its local atoms
    for (int i = 0; i < nlocal; i++) {
        int64_t j   = tag[i] - 1;
        double *x_i = x_local[i];
        double *x_g = &x[j * 3];
        EXPECT_DOUBLE_EQ(x_g[0], x_i[0]);
        EXPECT_DOUBLE_EQ(x_g[1], x_i[1]);
        EXPECT_DOUBLE_EQ(x_g[2], x_i[2]);
    }

    delete[] x;
}

TEST_F(MPITest, scatter)
{
    int *p_nlocal    = (int *)lammps_extract_global(lmp, "nlocal");
    int nlocal       = *p_nlocal;
    double *x_orig   = new double[3 * nlocal];
    double **x_local = (double **)lammps_extract_atom(lmp, "x");

    // make copy of original local x vector
    for (int i = 0; i < nlocal; i++) {
        int j         = 3 * i;
        x_orig[j]     = x_local[i][0];
        x_orig[j + 1] = x_local[i][1];
        x_orig[j + 2] = x_local[i][2];
    }

    // get the entire x on all procs
    int64_t natoms = (int64_t)lammps_get_natoms(lmp);
    double *x      = new double[natoms * 3];
    lammps_gather(lmp, (char *)"x", 1, 3, x);

    // shift all coordinates by 0.001
    const double delta = 0.001;
    for (int64_t i = 0; i < 3 * natoms; i++)
        x[i] += delta;

    // update positions of all atoms
    lammps_scatter(lmp, (char *)"x", 1, 3, x);
    delete[] x;
    x = nullptr;

    // get new nlocal and x_local
    p_nlocal = (int *)lammps_extract_global(lmp, "nlocal");
    nlocal   = *p_nlocal;
    x_local  = (double **)lammps_extract_atom(lmp, "x");

    ASSERT_EQ(nlocal, 8);

    // each proc checks its local atoms for shift
    for (int i = 0; i < nlocal; i++) {
        double *x_a = x_local[i];
        double *x_b = &x_orig[i * 3];
        EXPECT_DOUBLE_EQ(x_a[0], x_b[0] + delta);
        EXPECT_DOUBLE_EQ(x_a[1], x_b[1] + delta);
        EXPECT_DOUBLE_EQ(x_a[2], x_b[2] + delta);
    }

    delete[] x_orig;
}
#endif
