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

// Run the examples/ELASTIC workflow, which computes the zero temperature
// elastic stiffness tensor of Stillinger-Weber silicon from the stress
// response to finite box deformations.  The example is small enough (192
// atoms, minimizations capped at 100 iterations) to run unmodified, and
// its results can be checked against the analytical values for this
// potential quoted in the input: C11 = 151.4, C12 = 76.4, C44 = 56.4 GPa
// (E. R. Cowley, 1988).

#include "example_tests.h"

namespace LAMMPS_NS {

class ElasticExamplesTest : public ExampleTest {
protected:
    void SetUp() override
    {
        testbinary = "ElasticExamplesTest";
        LAMMPSTest::SetUp();
    }
};

TEST_F(ElasticExamplesTest, silicon_sw)
{
    REQUIRE_STYLES({"pair", "sw"});

    // in.elastic pulls in init.mod, potential.mod, and displace.mod with
    // relative paths and reads Si.sw from the same folder
    copy_from_examples("in.elastic");
    copy_from_examples("init.mod");
    copy_from_examples("potential.mod");
    copy_from_examples("displace.mod");
    copy_from_examples("Si.sw");

    BEGIN_HIDE_OUTPUT();
    command("include in.elastic");
    END_HIDE_OUTPUT();

    // cubic averages against the analytical values
    EXPECT_NEAR(get_variable_value("C11cubic"), 151.4, 1.5);
    EXPECT_NEAR(get_variable_value("C12cubic"), 76.4, 1.5);
    EXPECT_NEAR(get_variable_value("C44cubic"), 56.4, 1.5);
    EXPECT_NEAR(get_variable_value("bulkmodulus"), 101.4, 1.5);
    EXPECT_NEAR(get_variable_value("poissonratio"), 0.335, 0.01);

    // cubic symmetry: no coupling of shear and axial deformations
    EXPECT_NEAR(get_variable_value("C14all"), 0.0, 0.5);
    EXPECT_NEAR(get_variable_value("C25all"), 0.0, 0.5);
    EXPECT_NEAR(get_variable_value("C36all"), 0.0, 0.5);
    EXPECT_NEAR(get_variable_value("C45all"), 0.0, 0.5);

    delete_file("restart.equil");
    delete_file("in.elastic");
    delete_file("init.mod");
    delete_file("potential.mod");
    delete_file("displace.mod");
    delete_file("Si.sw");
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
