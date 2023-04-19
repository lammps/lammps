/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/core.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class FixMCBondTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "FixMCBondTest";
        LAMMPSTest::SetUp();
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in_double_bond.fourmol\"");
            command("group allwater molecule 3:6");
            command("region half block 0.0 INF INF INF INF INF");
            END_HIDE_OUTPUT();
        }
    }

    int get_num_bonds() { return (int)lammps_get_thermo(lmp, "bonds"); }

    // double get_scalar(const char *id)
    // {
    //     return *(double *)lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);
    // }

    // double *get_vector(const char *id)
    // {
    //     return (double *)lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    // }

    // double **get_array(const char *id)
    // {
    //     return (double **)lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY);
    // }
};

TEST_F(FixMCBondTest, BreakBond)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();
    int num_bonds_before = this->get_num_bonds();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");

    command("fix 1 bond/break 1 1 0.01");
    command("run 5");
    END_HIDE_OUTPUT();

    EXPECT_EQ(num_bonds_before, 27);
    EXPECT_EQ(this->get_num_bonds(), 27 - 3);
    // TODO: check for consistency of the inner structure of LAMMPS (i.e., special lists, angles
    // etc.)

    TEST_FAILURE(".*ERROR: Illegal fix bond/break command.*", command("fix 2 bond/break 1"););
    TEST_FAILURE(".*ERROR: Reuse of fix ID '1'.*", command("fix 1 bond/create 1"););
}

TEST_F(FixMCBondTest, CreateBond)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    int num_bonds_before = this->get_num_bonds();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");

    command("fix 1 bond/create 1 1 2 3.0 9");
    command("run 5");
    END_HIDE_OUTPUT();

    EXPECT_EQ(num_bonds_before, 27);
    EXPECT_GT(this->get_num_bonds(), num_bonds_before);
    // TODO: check for consistency of the inner structure of LAMMPS
    // (i.e., special lists, angles etc.)
}

TEST_F(FixMCBondTest, BreakSelf)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    int num_bonds_before = this->get_num_bonds();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");

    command("fix 1 bond/break/self 1 7");
    command("run 5");
    END_HIDE_OUTPUT();

    EXPECT_EQ(num_bonds_before, 27);
    EXPECT_EQ(this->get_num_bonds(), 25);
    // TODO: check for consistency of the inner structure of LAMMPS
    // (i.e., special lists, angles etc.)
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (LAMMPS_NS::platform::mpi_vendor() == "Open MPI" && !Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
