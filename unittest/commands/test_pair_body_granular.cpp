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

#include "lammps.h"

#include "force.h"
#include "info.h"
#include "library.h"
#include "pair.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

// pair style body/rounded/polyhedron/granular accepts only the GRANULAR
// sub-models that are well defined for the contacts of rounded polyhedra
// and for a single friction force per pair of bodies, see check_model()

class PairBodyGranularTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "PairBodyGranularTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite"};
        LAMMPSTest::SetUp();
        if (!Info::has_package("BODY") || !Info::has_package("GRANULAR")) GTEST_SKIP();

        BEGIN_HIDE_OUTPUT();
        command("units lj");
        command("atom_style body rounded/polyhedron 1 8");
        command("region box block 0 6 0 6 0 6");
        command("create_box 2 box");
        command("pair_style body/rounded/polyhedron/granular 0.01 0.2");
        END_HIDE_OUTPUT();
    }
};

TEST_F(PairBodyGranularTest, supported)
{
    BEGIN_HIDE_OUTPUT();
    command("pair_coeff 1 1 hooke 100 2 damping velocity tangential linear_history 30 0 0.2 "
            "cohesion 2");
    command("pair_coeff 1 2 hertz 80 2 damping mass_velocity tangential linear_nohistory 0 0.2");
    command("pair_coeff 2 2 hertz/material 1e4 0.5 0.3 tangential linear_history 300 0.5 0.2");
    END_HIDE_OUTPUT();
    ASSERT_NE(lmp->force->pair, nullptr);
    ASSERT_EQ(lmp->force->pair->setflag[1][1], 1);
    ASSERT_EQ(lmp->force->pair->setflag[1][2], 1);
    ASSERT_EQ(lmp->force->pair->setflag[2][2], 1);
}

TEST_F(PairBodyGranularTest, normal)
{
    const std::string tangential = " tangential linear_history 30 0 0.2";

    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the normal model dmt.*",
                 command("pair_coeff * * dmt 1e4 2 0.3 0.1" + tangential););
    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the normal model jkr.*",
                 command("pair_coeff * * jkr 1e4 2 0.3 0.1" + tangential););
    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the normal model mdr.*",
                 command("pair_coeff * * mdr 5e6 0.4 1.9e5 2.0 0.5 0.5 damping mdr 1" +
                         tangential););
}

TEST_F(PairBodyGranularTest, damping)
{
    const std::string normal     = "pair_coeff * * hertz/material 1e4 0.5 0.3";
    const std::string tangential = " tangential linear_history 300 0.5 0.2";

    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the damping model viscoelastic.*",
                 command(normal + " damping viscoelastic" + tangential););
    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the damping model tsuji.*",
                 command(normal + " damping tsuji" + tangential););
    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the damping model coeff_restitution.*",
                 command(normal + " damping coeff_restitution" + tangential););
}

TEST_F(PairBodyGranularTest, tangential)
{
    const std::string normal = "pair_coeff * * hertz/material 1e4 0.5 0.3 tangential ";

    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the tangential model mindlin.*",
                 command(normal + "mindlin 300 0.5 0.2"););
    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the tangential model mindlin/force.*",
                 command(normal + "mindlin/force 300 0.5 0.2"););
    TEST_FAILURE(".*ERROR: Pair style body/rounded/polyhedron/granular does not support "
                 "the tangential model mindlin_rescale.*",
                 command(normal + "mindlin_rescale 300 0.5 0.2"););
}

TEST_F(PairBodyGranularTest, keywords)
{
    // rolling, twisting, heat, and limit_damping are not keywords of this pair style

    const std::string model = "pair_coeff * * hooke 100 2 tangential linear_history 30 0 0.2 ";

    TEST_FAILURE(".*ERROR: Unknown pair_coeff keyword rolling for pair style "
                 "body/rounded/polyhedron/granular.*",
                 command(model + "rolling sds 100 0 0.1"););
    TEST_FAILURE(".*ERROR: Unknown pair_coeff keyword twisting for pair style "
                 "body/rounded/polyhedron/granular.*",
                 command(model + "twisting marshall"););
    TEST_FAILURE(".*ERROR: Unknown pair_coeff keyword heat for pair style "
                 "body/rounded/polyhedron/granular.*",
                 command(model + "heat area 1.0"););
    TEST_FAILURE(".*ERROR: Unknown pair_coeff keyword limit_damping for pair style "
                 "body/rounded/polyhedron/granular.*",
                 command(model + "limit_damping"););
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
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
