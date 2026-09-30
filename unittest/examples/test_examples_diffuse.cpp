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

// Run abbreviated versions of the examples/DIFFUSE inputs, which compute
// the diffusion coefficient of a 2d LJ liquid from the mean-squared
// displacement (compute msd) and from integrating the velocity auto-
// correlation function (compute vacf).  The truncated runs cannot deliver
// converged values, but a liquid must yield a positive estimate from
// either method.

#include "example_tests.h"

namespace LAMMPS_NS {

class DiffuseExamplesTest : public ExampleTest {
protected:
    void SetUp() override
    {
        testbinary = "DiffuseExamplesTest";
        LAMMPSTest::SetUp();
        // abbreviate: 200 atoms instead of 3200, short equilibration and
        // data gathering runs
        preset("x", "10");
        preset("y", "10");
        preset("nequil", "200");
        preset("nprod", "500");
    }
};

TEST_F(DiffuseExamplesTest, msd)
{
    run_input("in.msd.2d");
    ASSERT_EQ(lmp->update->ntimestep, 500);

    // the mean-squared displacement of a liquid grows with time, so both
    // estimates of the diffusion coefficient must be positive: v_twopoint
    // is MSD/(4t) at the final step, v_fitslope a fit to MSD vs. time
    const double msd      = get_variable_value("twopoint");
    const double fitslope = get_variable_value("fitslope");
    EXPECT_TRUE(std::isfinite(msd));
    EXPECT_TRUE(std::isfinite(fitslope));
    EXPECT_GT(msd, 0.0);
    EXPECT_GT(fitslope, 0.0);
}

TEST_F(DiffuseExamplesTest, vacf)
{
    run_input("in.vacf.2d");
    ASSERT_EQ(lmp->update->ntimestep, 500);

    // v_vacf integrates the velocity autocorrelation function; over the
    // short run the positive short-time correlation dominates
    const double vacf = get_variable_value("vacf");
    EXPECT_TRUE(std::isfinite(vacf));
    EXPECT_GT(vacf, 0.0);
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
