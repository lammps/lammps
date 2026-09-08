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

// Run abbreviated versions of the examples/ELASTIC_T workflows, which
// compute the elastic stiffness tensor at finite temperature: DEFORMATION
// measures the change of the average stress tensor under finite box
// deformations, BORN_MATRIX averages the Born matrix and the stress
// fluctuations of a single trajectory (compute born/matrix, analytically
// or via numerical differentiation).  The truncated sampling cannot give
// converged elastic constants, so the checks are limited to the workflows
// completing with plausible values.

#include "example_tests.h"

namespace LAMMPS_NS {

class ElasticTExamplesTest : public ExampleTest {
protected:
    void SetUp() override
    {
        testbinary = "ElasticTExamplesTest";
        LAMMPSTest::SetUp();
    }
};

TEST_F(ElasticTExamplesTest, deformation_silicon)
{
    REQUIRE_STYLES({"pair", "sw"});

    copy_from_examples("DEFORMATION/Silicon/in.elastic");
    copy_from_examples("DEFORMATION/Silicon/init.mod");
    copy_from_examples("DEFORMATION/Silicon/potential.mod");
    copy_from_examples("DEFORMATION/Silicon/displace.mod");
    copy_from_examples("DEFORMATION/Silicon/Si.sw");

    // abbreviate: 10 stress samples per average instead of 100, which
    // shortens all run segments accordingly (they are multiples of the
    // averaging window)
    preset("nevery", "2");
    preset("nrepeat", "5");

    BEGIN_HIDE_OUTPUT();
    command("include in.elastic");
    END_HIDE_OUTPUT();

    // at 2000 K with this little sampling only the sign and magnitude of
    // the diagonal averages are robust
    const double c11  = get_variable_value("C11cubic");
    const double c44  = get_variable_value("C44cubic");
    const double bulk = get_variable_value("bulkmodulus");
    EXPECT_TRUE(std::isfinite(c11));
    EXPECT_TRUE(std::isfinite(c44));
    EXPECT_TRUE(std::isfinite(bulk));
    EXPECT_GT(c11, 0.0);
    EXPECT_LT(c11, 500.0);
    EXPECT_GT(bulk, 0.0);

    delete_file("restart.equil");
    delete_file("in.elastic");
    delete_file("init.mod");
    delete_file("potential.mod");
    delete_file("displace.mod");
    delete_file("Si.sw");
}

TEST_F(ElasticTExamplesTest, born_matrix_argon)
{
    REQUIRE_STYLES({"compute", "born/matrix"});

    // abbreviate equilibration and sampling
    preset("nequil", "200");
    preset("nsteps", "300");
    preset("nthermo", "100");

    // the analytical and the numerical differentiation Born matrix run the
    // same trajectory (identical dynamics and seeds), so their elastic
    // constants must agree to within the numerical differentiation error
    run_input("BORN_MATRIX/Argon/Analytical/in.ljcov");
    const double c11_ana = get_variable_value("aC11");
    const double c12_ana = get_variable_value("aC12");
    const double c44_ana = get_variable_value("aC44");
    EXPECT_TRUE(std::isfinite(c11_ana));
    EXPECT_TRUE(std::isfinite(c12_ana));
    EXPECT_TRUE(std::isfinite(c44_ana));

    run_input("BORN_MATRIX/Argon/Numdiff/in.ljcov");
    const double c11_num = get_variable_value("aC11");
    const double c12_num = get_variable_value("aC12");
    const double c44_num = get_variable_value("aC44");

    EXPECT_NEAR(c11_ana, c11_num, 0.01 * std::fabs(c11_ana) + 0.1);
    EXPECT_NEAR(c12_ana, c12_num, 0.01 * std::fabs(c12_ana) + 0.1);
    EXPECT_NEAR(c44_ana, c44_num, 0.01 * std::fabs(c44_ana) + 0.1);

    delete_file("born.out");
    delete_file("vir.out");
}

TEST_F(ElasticTExamplesTest, born_matrix_silicon)
{
    REQUIRE_STYLES({"pair", "sw"}, {"compute", "born/matrix"});

    copy_from_examples("BORN_MATRIX/Silicon/in.elastic");
    copy_from_examples("BORN_MATRIX/Silicon/init.in");
    copy_from_examples("BORN_MATRIX/Silicon/potential.in");
    copy_from_examples("BORN_MATRIX/Silicon/output.in");
    copy_from_examples("BORN_MATRIX/Silicon/final_output.in");
    copy_from_examples("BORN_MATRIX/Silicon/Si.sw");

    // abbreviate: the averaging window drives all run lengths
    preset("nthermo", "50");
    preset("neveryborn", "10");
    preset("logsuffix", "test");

    BEGIN_HIDE_OUTPUT();
    command("include in.elastic");
    END_HIDE_OUTPUT();

    EXPECT_TRUE(std::isfinite(get_variable_value("C11cubic")));
    EXPECT_TRUE(std::isfinite(get_variable_value("C12cubic")));
    EXPECT_TRUE(std::isfinite(get_variable_value("C44cubic")));

    delete_file("log.elastic.test");
    delete_file("in.elastic");
    delete_file("init.in");
    delete_file("potential.in");
    delete_file("output.in");
    delete_file("final_output.in");
    delete_file("Si.sw");
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
