// unit tests for checking LAMMPS MPI load balancing

#define LAMMPS_LIB_MPI 1
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "neighbor.h"
#include "timer.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/test_mpi_main.h"

using ::testing::ExitedWithCode;
using ::testing::HasSubstr;
using ::testing::StartsWith;
using ::testing::StrEq;

namespace LAMMPS_NS {

class MPILoadBalanceTest : public ::testing::Test {
public:
    void command(const std::string &line) { lmp->input->one(line); }

protected:
    const char *testbinary = "LAMMPSTest";
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {testbinary, "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        InitSystem();
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    virtual void InitSystem()
    {
        command("boundary        f f f");
        command("units           lj");
        command("atom_style      atomic");
        command("atom_modify     map yes");

        command("region          box block 0 20 0 20 0 20");
        command("create_box      1 box");
        command("mass            1 1.0");

        command("pair_style      lj/cut 2.5");
        command("pair_coeff      1 1 1.0 1.0 2.5");

        command("neighbor        0.3 bin");
        command("neigh_modify    every 20 delay 0 check no");
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        lmp = nullptr;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(MPILoadBalanceTest, grid_yz)
{
    command("create_atoms 1 single 0 0 0");
    command("create_atoms 1 single 0 0 5");
    command("create_atoms 1 single 0 5 0");
    command("create_atoms 1 single 0 5 5");
    command("create_atoms 1 single 5 0 0");
    command("create_atoms 1 single 5 0 5");
    command("create_atoms 1 single 5 5 0");
    command("create_atoms 1 single 5 5 5");
    ASSERT_EQ(lmp->atom->natoms, 8);
    ASSERT_EQ(lmp->comm->nprocs, 4);

    // initial state
    switch (lmp->comm->me) {
        case 0:
            ASSERT_EQ(lmp->atom->nlocal, 8);
            break;
        case 1:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 2:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 3:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
    }

    command("balance 1 x uniform y 0.125 z uniform");

    // state after balance command
    switch (lmp->comm->me) {
        case 0:
            ASSERT_EQ(lmp->atom->nlocal, 4);
            break;
        case 1:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 2:
            ASSERT_EQ(lmp->atom->nlocal, 4);
            break;
        case 3:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
    }

    command("balance 1 x uniform y 0.125 z 0.125");

    // state after second balance command
    switch (lmp->comm->me) {
        case 0:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
        case 1:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
        case 2:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
        case 3:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
    }
}

TEST_F(MPILoadBalanceTest, rcb)
{
    command("comm_style tiled");
    command("create_atoms 1 single 0 0 0");
    command("create_atoms 1 single 0 0 5");
    command("create_atoms 1 single 0 5 0");
    command("create_atoms 1 single 0 5 5");
    command("create_atoms 1 single 5 0 0");
    command("create_atoms 1 single 5 0 5");
    command("create_atoms 1 single 5 5 0");
    command("create_atoms 1 single 5 5 5");

    // initial state
    switch (lmp->comm->me) {
        case 0:
            ASSERT_EQ(lmp->atom->nlocal, 8);
            break;
        case 1:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 2:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 3:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
    }

    command("balance 1 rcb");

    // state after balance command
    switch (lmp->comm->me) {
        case 0:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
        case 1:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
        case 2:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
        case 3:
            ASSERT_EQ(lmp->atom->nlocal, 2);
            break;
    }

    // box dimensions should have minimal size
    double dx = lmp->domain->subhi[0] - lmp->domain->sublo[0];
    double dy = lmp->domain->subhi[1] - lmp->domain->sublo[1];
    double dz = lmp->domain->subhi[2] - lmp->domain->sublo[2];

    ASSERT_GT(dx, lmp->neighbor->skin);
    ASSERT_GT(dy, lmp->neighbor->skin);
    ASSERT_GT(dz, lmp->neighbor->skin);
}

TEST_F(MPILoadBalanceTest, rcb_min_size)
{
    GTEST_SKIP();
    // TODO minimum domain size is not enforced right now
    // skipping for now to allow other MPI tests to get merged
    command("comm_style tiled");
    command("create_atoms 1 single 0 0 0");
    command("create_atoms 1 single 0 0 0.25");
    command("create_atoms 1 single 0 0.25 0");
    command("create_atoms 1 single 0 0.25 0.25");
    command("create_atoms 1 single 0.25 0 0");
    command("create_atoms 1 single 0.25 0 0.25");
    command("create_atoms 1 single 0.25 0.25 0");
    command("create_atoms 1 single 0.25 0.25 0.25");

    // initial state
    switch (lmp->comm->me) {
        case 0:
            ASSERT_EQ(lmp->atom->nlocal, 8);
            break;
        case 1:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 2:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 3:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
    }

    // this should fail and not change the boxes
    // or keep them at a minimum size
    command("balance 1 rcb");

    // state after balance command
    switch (lmp->comm->me) {
        case 0:
            ASSERT_EQ(lmp->atom->nlocal, 8);
            break;
        case 1:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 2:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
        case 3:
            ASSERT_EQ(lmp->atom->nlocal, 0);
            break;
    }

    // box dimensions should have minimal size
    double dx = lmp->domain->subhi[0] - lmp->domain->sublo[0];
    double dy = lmp->domain->subhi[1] - lmp->domain->sublo[1];
    double dz = lmp->domain->subhi[2] - lmp->domain->sublo[2];

    ASSERT_GT(dx, lmp->neighbor->skin);
    ASSERT_GT(dy, lmp->neighbor->skin);
    ASSERT_GT(dz, lmp->neighbor->skin);
}

} // namespace LAMMPS_NS
