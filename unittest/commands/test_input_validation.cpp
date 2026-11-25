/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Regression tests for input values that used to be silently accepted
// and later caused reads of invalid or uninitialized data.

#include "lammps.h"

#include "info.h"
#include "library.h"
#include "utils.h"

#include "../testing/core.h"
#include "../testing/utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <fstream>
#include <vector>

using namespace LAMMPS_NS;

bool verbose = false;

class InputValidationTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "InputValidationTest";
        LAMMPSTest::SetUp();
    }

    void atomic_box()
    {
        BEGIN_HIDE_OUTPUT();
        command("units lj");
        command("region box block 0 2 0 2 0 2");
        command("create_box 1 box");
        command("create_atoms 1 single 1.0 1.0 1.0");
        command("mass 1 1.0");
        END_HIDE_OUTPUT();
    }
};

TEST_F(InputValidationTest, improper_cvff_multiplicity)
{
    if (!Info::has_package("MOLECULE")) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("units real");
    command("atom_style full");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box improper/types 1 extra/improper/per/atom 1");
    command("improper_style cvff");
    END_HIDE_OUTPUT();

    // the harmonic lookup tables only cover multiplicities 0 to 6

    TEST_FAILURE(".*ERROR: Improper cvff multiplicity 7 must be between 0 and 6.*",
                 command("improper_coeff 1 80.0 1 7"););
    TEST_FAILURE(".*ERROR: Improper cvff multiplicity -1 must be between 0 and 6.*",
                 command("improper_coeff 1 80.0 1 -1"););
    BEGIN_HIDE_OUTPUT();
    command("improper_coeff 1 80.0 1 2");
    END_HIDE_OUTPUT();
}

TEST_F(InputValidationTest, fix_vector_nmax)
{
    atomic_box();
    BEGIN_HIDE_OUTPUT();
    command("compute t all temp");
    END_HIDE_OUTPUT();

    // nmax values beyond 32-bit int would truncate the allocation

    TEST_FAILURE(".*ERROR: Invalid nmax value.*",
                 command("fix v1 all vector 1 c_t nmax 3000000000"););
    TEST_FAILURE(".*ERROR: Invalid nmax value.*", command("fix v2 all vector 1 c_t nmax 0"););
    BEGIN_HIDE_OUTPUT();
    command("fix v3 all vector 1 c_t nmax 1000");
    END_HIDE_OUTPUT();
}

TEST_F(InputValidationTest, fix_ave_grid_indexed_variable)
{
    atomic_box();
    BEGIN_HIDE_OUTPUT();
    command("variable f atom x+y");
    END_HIDE_OUTPUT();

    // atom-style variables produce a vector: a bracketed reference has
    // no array to read and used to access uninitialized storage

    TEST_FAILURE(".*ERROR: Fix ave/grid variable f cannot be indexed.*",
                 command("fix a1 all ave/grid 1 1 1 2 2 2 v_f[2]"););
    BEGIN_HIDE_OUTPUT();
    command("fix a2 all ave/grid 1 1 1 2 2 2 v_f");
    END_HIDE_OUTPUT();
}

TEST_F(InputValidationTest, fix_store_state_no_atoms)
{
    // atom-style variable input with zero atoms on a proc used to
    // dereference the (then null) per-atom storage array

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("region box block 0 2 0 2 0 2");
    command("create_box 1 box");
    command("mass 1 1.0");
    command("variable f atom x+y");
    command("fix s all store/state 1 v_f");
    command("run 1 post no");
    END_HIDE_OUTPUT();
    SUCCEED();
}

TEST_F(InputValidationTest, pair_tersoff_empty_potential)
{
    if (!Info::has_package("MANYBODY")) GTEST_SKIP();

    // an empty or comment-only potential file used to hang forever

    {
        std::ofstream out("empty_test.tersoff");
        out << "# DATE: 2007-06-11\n\n";
    }
    atomic_box();
    BEGIN_HIDE_OUTPUT();
    command("pair_style tersoff");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*Potential file is missing an entry for.*",
                 command("pair_coeff * * empty_test.tersoff X"););
    delete_file("empty_test.tersoff");
}

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
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
