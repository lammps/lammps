// unit tests for checking LAMMPS MPI load balancing

#define LAMMPS_LIB_MPI 1
#include "lammps.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "input.h"
#include "timer.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/test_mpi_main.h"

using ::testing::ExitedWithCode;
using ::testing::HasSubstr;
using ::testing::StartsWith;
using ::testing::StrEq;

namespace LAMMPS_NS
{

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
        command("units           lj");
        command("atom_style      atomic");
        command("atom_modify     map yes");

        command("region          box block 0 2 0 2 0 2");
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
    command("create_atoms 1 single 0 0 0.5");
    command("create_atoms 1 single 0 0.5 0");
    command("create_atoms 1 single 0 0.5 0.5");
    command("create_atoms 1 single 0.5 0 0");
    command("create_atoms 1 single 0.5 0 0.5");
    command("create_atoms 1 single 0.5 0.5 0");
    command("create_atoms 1 single 0.5 0.5 0.5");
    ASSERT_EQ(lmp->atom->natoms, 8);
    ASSERT_EQ(lmp->comm->nprocs, 4);

    switch(lmp->comm->me) {
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

    command("balance 1 x uniform y 0.25 z uniform");

    switch(lmp->comm->me) {
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

    command("balance 1 x uniform y 0.25 z 0.25");

    switch(lmp->comm->me) {
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

}
