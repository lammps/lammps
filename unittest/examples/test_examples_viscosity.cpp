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

// Run abbreviated versions of the examples/VISCOSITY inputs, which compute
// the shear viscosity of an LJ liquid (or SPC/E water for the cosine
// acceleration method) with different methods.  The truncated runs cannot
// give converged values: where the estimate is positive by construction
// (Muller-Plathe momentum flux, NEMD shear response, Einstein sum of
// squares, driven cosine profile) the checks require that, otherwise only
// finiteness.

#include "example_tests.h"

namespace LAMMPS_NS {

class ViscosityExamplesTest : public ExampleTest {
protected:
    void SetUp() override
    {
        testbinary = "ViscosityExamplesTest";
        LAMMPSTest::SetUp();
        // abbreviate: 200 atoms instead of 800 for the 2d LJ inputs
        preset("x", "10");
        preset("y", "10");
    }
};

TEST_F(ViscosityExamplesTest, muller_plathe)
{
    REQUIRE_STYLES({"fix", "viscosity"});
    preset("nequil", "1000");
    preset("nsteady", "1000");
    preset("nprod", "2000");
    run_input("in.mp.2d");

    // the momentum flux imposed by fix viscosity and the resulting
    // velocity profile difference have opposite signs by construction
    const double visc = get_variable_value("visc");
    EXPECT_TRUE(std::isfinite(visc));
    EXPECT_GT(visc, 0.0);

    auto profile = last_vector_block("profile.mp.2d", 3);
    ASSERT_EQ(profile.size(), 20);
    delete_file("profile.mp.2d");
}

TEST_F(ViscosityExamplesTest, nemd_sllod)
{
    preset("nequil", "1000");
    preset("nsteady", "5000");
    preset("nprod", "1000");
    preset("crepeat", "50");
    preset("cfreq", "1000");
    run_input("in.nemd.2d");

    // under imposed shear the average of -pxy/srate must be positive
    BEGIN_HIDE_OUTPUT();
    command("variable vave_check equal f_vave");
    END_HIDE_OUTPUT();
    const double vave = get_variable_value("vave_check");
    EXPECT_TRUE(std::isfinite(vave));
    EXPECT_GT(vave, 0.0);

    auto profile = last_vector_block("profile.nemd.2d", 3);
    ASSERT_EQ(profile.size(), 20);
    delete_file("profile.nemd.2d");
}

TEST_F(ViscosityExamplesTest, einstein)
{
    preset("p", "40");
    preset("s", "5");
    preset("nequil", "1000");
    preset("nprod", "1000");
    run_input("in.einstein.2d");

    // sums of squared stress displacements scaled by a positive factor
    const double eta = get_variable_value("eta");
    EXPECT_TRUE(std::isfinite(eta));
    EXPECT_GT(eta, 0.0);
    delete_file("profile.einstein.2d");
}

TEST_F(ViscosityExamplesTest, green_kubo)
{
    preset("p", "40");
    preset("s", "5");
    preset("nequil", "1000");
    preset("nprod", "1000");
    run_input("in.gk.2d");

    // the Green-Kubo integral fluctuates in sign for short sampling
    const double eta = get_variable_value("eta");
    EXPECT_TRUE(std::isfinite(eta));
    ASSERT_FILE_EXISTS("profile.gk.2d");
    delete_file("profile.gk.2d");
}

TEST_F(ViscosityExamplesTest, shearing_wall)
{
    preset("nequil", "2000");
    preset("nprod", "1000");
    preset("wrepeat", "50");
    preset("wfreq", "1000");
    run_input("in.wall.2d");

    // the wall driven shear develops slowly, so the short time running
    // average can still have either sign
    BEGIN_HIDE_OUTPUT();
    command("variable vave_check equal f_vave");
    END_HIDE_OUTPUT();
    EXPECT_TRUE(std::isfinite(get_variable_value("vave_check")));

    auto profile = last_vector_block("profile.wall.2d", 3);
    ASSERT_EQ(profile.size(), 20);
    delete_file("profile.wall.2d");
}

TEST_F(ViscosityExamplesTest, cosine_acceleration)
{
    REQUIRE_STYLES({"fix", "accelerate/cos"}, {"compute", "viscosity/cos"}, {"fix", "shake"},
                   {"kspace", "pppm"}, {"bond", "harmonic"});
    copy_from_examples("data.cos.1000SPCE");
    preset("nprod", "250");
    run_input("in.cos.1000SPCE");

    // the cosine acceleration drives a velocity profile with positive
    // amplitude, so the reciprocal viscosity estimate must be positive
    const double invvis = get_variable_value("invVis");
    EXPECT_TRUE(std::isfinite(invvis));
    EXPECT_GT(invvis, 0.0);

    delete_file("data.cos.1000SPCE");
    delete_file("dump.lammpstrj");
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
