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

// Run abbreviated versions of the examples/HEAT inputs, which establish a
// thermal gradient in an LJ fluid or SPC/E water with the HEX or eHEX
// heat exchange algorithm of fix ehex (Wirnsberger, Frenkel, Dellago,
// JCP 143, 124104 (2015)).  The truncated runs only exercise the workflow:
// the checks are that the run completes with finite total energy and that
// the temperature profile is written with plausible values.

#include "example_tests.h"

namespace LAMMPS_NS {

class HeatExamplesTest : public ExampleTest {
protected:
    void check_profile(const std::string &filename)
    {
        // temperature profile from fix ave/chunk: 60 bins with positive,
        // finite temperatures in the fourth column
        ASSERT_FILE_EXISTS(filename);
        auto profile = last_vector_block(filename, 3);
        ASSERT_EQ(profile.size(), 60);
        for (const auto &row : profile) {
            ASSERT_EQ(row.size(), 4);
            EXPECT_TRUE(std::isfinite(row[3]));
            EXPECT_GE(row[3], 0.0);
        }
    }

    void run_lj(const std::string &variant)
    {
        REQUIRE_STYLES({"pair", "lj/sf"}, {"fix", "ehex"});
        copy_from_examples("data.lj");
        // abbreviate: ~1100 steps of production instead of ~714000
        preset("tprod", "7.7");
        run_input("in.lj." + variant);

        EXPECT_GT(lmp->update->ntimestep, 1000);
        EXPECT_TRUE(std::isfinite(thermo_value("etotal")));
        check_profile("out.Tlj_" + variant);

        delete_file("data.lj");
        delete_file("log.lj_" + variant);
        delete_file("out.Tlj_" + variant);
        delete_file("out.Elj_" + variant);
    }

    void run_spce(const std::string &variant)
    {
        REQUIRE_STYLES({"kspace", "ewald"}, {"bond", "harmonic"}, {"fix", "rattle"},
                       {"fix", "ehex"});
        copy_from_examples("data.spce");
        // abbreviate: 100 steps of production instead of ~333000
        preset("tprod", "300");
        run_input("in.spce." + variant);

        EXPECT_EQ(lmp->update->ntimestep, 100);
        EXPECT_TRUE(std::isfinite(thermo_value("etotal")));
        check_profile("out.Tspce_" + variant);

        delete_file("data.spce");
        delete_file("log.spce_" + variant);
        delete_file("out.Tspce_" + variant);
        delete_file("out.Espce_" + variant);
    }

    void SetUp() override
    {
        testbinary = "HeatExamplesTest";
        LAMMPSTest::SetUp();
    }
};

TEST_F(HeatExamplesTest, lj_hex)
{
    run_lj("hex");
}

TEST_F(HeatExamplesTest, lj_ehex)
{
    run_lj("ehex");
}

TEST_F(HeatExamplesTest, spce_hex)
{
    run_spce("hex");
}

TEST_F(HeatExamplesTest, spce_ehex)
{
    run_spce("ehex");
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
