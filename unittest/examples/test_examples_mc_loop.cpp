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

// Run an abbreviated version of the examples/MC-LOOP input, which uses
// LAMMPS as the energy evaluation engine of a Metropolis Monte Carlo
// relaxation loop implemented with input script commands (label/jump,
// variables, run 1 pre no post no).  The loop count is reduced from 3000
// to 30 moves; a loop-style variable definition in the input is skipped
// when the preset variable of the same name already exists.

#include "example_tests.h"

namespace LAMMPS_NS {

class MCLoopExamplesTest : public ExampleTest {
protected:
    void SetUp() override
    {
        testbinary = "MCLoopExamplesTest";
        LAMMPSTest::SetUp();
        HIDE_OUTPUT([&] {
            command("variable iter loop 30");
        });
    }
};

TEST_F(MCLoopExamplesTest, metropolis_relaxation)
{
    run_input("in.mc");

    const double emin    = get_variable_value("emin");
    const double estart  = get_variable_value("estart");
    const double efinal  = get_variable_value("e");
    const double naccept = get_variable_value("naccept");

    // the perturbed lattice starts well above the perfect lattice energy
    EXPECT_GT(estart, emin);
    // relaxation must not go below the perfect lattice minimum and moves
    // are only accepted downhill or with a tiny Boltzmann factor uphill
    EXPECT_GE(efinal, emin - 1.0e-10 * std::fabs(emin));
    EXPECT_LE(efinal, estart + 1.0e-10 * std::fabs(estart));
    EXPECT_GE(naccept, 0.0);
    EXPECT_LE(naccept, 30.0);
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
