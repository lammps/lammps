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

// Run abbreviated versions of the examples/KAPPA inputs, which compute the
// thermal conductivity of an LJ liquid with 5 different methods.  The
// truncated runs cannot give converged values, but for the thermal
// gradient methods the estimate must come out positive because the sign
// of the numerator (exchanged energy) and denominator (temperature
// difference) is fixed by construction; for the Green-Kubo method only
// finiteness can be required.

#include "example_tests.h"

namespace LAMMPS_NS {

class KappaExamplesTest : public ExampleTest {
protected:
    void SetUp() override
    {
        testbinary = "KappaExamplesTest";
        LAMMPSTest::SetUp();
        // abbreviate: 720 instead of 8000 atoms.  the z extent stays at 20
        // lattice cells since the heated/cooled regions and the reported
        // bin indices are tied to it.  the run lengths are preset per test,
        // they must stay commensurate with the Nfreq = 1000 averaging
        // windows of the inputs
        preset("x", "3");
        preset("y", "3");
    }

    void check_kappa(const std::string &profile, bool positive = true)
    {
        const double kappa = get_variable_value("kappa");
        EXPECT_TRUE(std::isfinite(kappa));
        if (positive) EXPECT_GT(kappa, 0.0);

        if (!profile.empty()) {
            // temperature profile with 20 bins from fix ave/chunk
            ASSERT_FILE_EXISTS(profile);
            auto block = last_vector_block(profile, 3);
            ASSERT_EQ(block.size(), 20);
            for (const auto &row : block) {
                ASSERT_EQ(row.size(), 4);
                EXPECT_TRUE(std::isfinite(row[3]));
                EXPECT_GE(row[3], 0.0);
            }
            delete_file(profile);
        }
    }
};

TEST_F(KappaExamplesTest, muller_plathe)
{
    // the thermo output of the second equilibration already reports the
    // temperature profile difference, so the equilibration segments must
    // end on the Nfreq = 1000 averaging grid
    preset("nequil", "1000");
    preset("nsteady", "1000");
    preset("nprod", "2000");
    run_input("in.mp");
    check_kappa("profile.mp");
}

TEST_F(KappaExamplesTest, fix_heat)
{
    preset("nequil", "500");
    preset("nsteady", "500");
    preset("nprod", "2000");
    run_input("in.heat");
    check_kappa("profile.heat");
}

TEST_F(KappaExamplesTest, fix_ehex)
{
    REQUIRE_STYLES({"fix", "ehex"});
    preset("nequil", "500");
    preset("nsteady", "500");
    preset("nprod", "2000");
    run_input("in.ehex");
    check_kappa("profile.ehex");
}

TEST_F(KappaExamplesTest, langevin)
{
    preset("nequil", "500");
    preset("nsteady", "500");
    preset("nprod", "2000");
    run_input("in.langevin");
    check_kappa("profile.langevin");
}

TEST_F(KappaExamplesTest, green_kubo_heatflux)
{
    // different structure: no thermal gradient, correlation windows instead
    preset("p", "20");
    preset("s", "5");
    preset("nequil", "200");
    preset("nprod", "1000");
    run_input("in.heatflux");
    // the Green-Kubo integral fluctuates in sign for short sampling
    check_kappa("", false);
    ASSERT_FILE_EXISTS("profile.heatflux");
    delete_file("profile.heatflux");
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
