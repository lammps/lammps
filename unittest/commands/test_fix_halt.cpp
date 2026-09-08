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
#include "lammps.h"
#include "library.h"
#include "update.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

using ::testing::ContainsRegex;

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// fix halt stops a run when a criterion is met, so it cannot be checked with
// the fixed step count drivers in unittest/force-styles.  the tests below run
// the same input with the plain and - through the LAMMPS_ACCELERATOR_ARGS
// setting of the FixHaltKokkos ctest entry - with the KOKKOS version of the
// fix and confirm that the run stops on the expected timestep.

class FixHaltTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "FixHaltTest";
        LAMMPSTest::SetUp();
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in.fourmol\"");
            command("fix move all nve");
            END_HIDE_OUTPUT();
        }
    }
};

// the longest bond of the fourmol system is about 1.565 angstrom and the pair
// and bond styles of the input template are "zero", so the geometry hardly
// changes during the short run below.  a limit of 1.5 angstrom is thus already
// exceeded at the first end-of-step check.

TEST_F(FixHaltTest, BondmaxTriggers)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("fix stop all halt 1 bondmax > 1.5");
    END_HIDE_OUTPUT();
    auto output = CAPTURE_OUTPUT([&] {
        command("run 10 post no");
    });

    EXPECT_EQ(lmp->update->ntimestep, 1);
    EXPECT_THAT(output, ContainsRegex("Fix halt condition for fix-id stop met on step 1"));
}

// with a limit well above the longest bond the run must not be stopped early

TEST_F(FixHaltTest, BondmaxDoesNotTrigger)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("fix stop all halt 1 bondmax > 2.0");
    command("run 10 post no");
    END_HIDE_OUTPUT();

    EXPECT_EQ(lmp->update->ntimestep, 10);
}

// an equal-style variable is only checked every "nevery" steps, so the run
// stops on step 6 and not on step 5, where the condition first becomes true

TEST_F(FixHaltTest, VariableHalt)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("variable steps equal step");
    command("fix stop all halt 2 v_steps >= 5");
    END_HIDE_OUTPUT();
    auto output = CAPTURE_OUTPUT([&] {
        command("run 20 post no");
    });

    EXPECT_EQ(lmp->update->ntimestep, 6);
    EXPECT_THAT(output, ContainsRegex("Fix halt condition for fix-id stop met on step 6"));
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") verbose = true;
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
