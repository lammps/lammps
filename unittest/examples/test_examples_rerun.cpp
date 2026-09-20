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

// Run the input scripts of the examples/rerun folder in the order the README
// there prescribes and check the claims it makes.  The box length is reduced
// from 20 to 6 lattice cells (864 instead of 32000 atoms): the workflow only
// needs to be seen at work, not converged.

#include "example_tests.h"

namespace LAMMPS_NS {

class RerunExamplesTest : public ExampleTest {
protected:
    void SetUp() override
    {
        testbinary = "RerunExamplesTest";
        LAMMPSTest::SetUp();
        // abbreviate: 6^3 fcc cells instead of 20^3.  the box must still be
        // at least twice the 5 sigma g(r) cutoff of in.rdf.rerun across
        preset("len", "6");
    }
};

TEST_F(RerunExamplesTest, first_rerun_read_dump)
{
    // writes the snapshots 0 to 1000 of an LJ melt to lj.dump
    run_input("in.first");
    ASSERT_FILE_EXISTS("lj.dump");
    ASSERT_EQ(lmp->update->ntimestep, 1000);
    const double pe_live   = thermo_value("pe");
    const double temp_live = thermo_value("temp");

    // re-read the final snapshot in the same instance: the differences to
    // the live run come only from the "%g" precision of the dump file
    BEGIN_HIDE_OUTPUT();
    command("undump 1");
    command("read_dump lj.dump 1000 x y z vx vy vz");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    EXPECT_NEAR(thermo_value("pe"), pe_live, 1.0e-5 * std::fabs(pe_live));
    EXPECT_NEAR(thermo_value("temp"), temp_live, 1.0e-5 * temp_live);

    // in.rerun recomputes thermo output on the snapshots 200 to 800
    run_input("in.rerun");
    ASSERT_EQ(lmp->update->ntimestep, 800);
    const double pe_rerun = thermo_value("pe");

    // the rerun and read_dump commands must produce the same state from the
    // same snapshot
    BEGIN_HIDE_OUTPUT();
    command("read_dump lj.dump 800 x y z vx vy vz");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    EXPECT_NEAR(thermo_value("pe"), pe_rerun, 1.0e-10 * std::fabs(pe_rerun));

    // in.read_dump visits the snapshots 200, 800, 600, and 400
    run_input("in.read_dump");
    ASSERT_EQ(lmp->update->ntimestep, 400);
    const double pe_read = thermo_value("pe");

    BEGIN_HIDE_OUTPUT();
    command("rerun lj.dump first 400 last 400 every 1 post no dump x y z vx vy vz");
    END_HIDE_OUTPUT();
    EXPECT_NEAR(thermo_value("pe"), pe_read, 1.0e-10 * std::fabs(pe_read));

    delete_file("lj.dump");
}

TEST_F(RerunExamplesTest, rdf_first_rerun)
{
    // writes its own lj.dump plus the time averaged RDF with 50 bins out to
    // a distance of 2.5 sigma to rdf.first
    run_input("in.rdf.first");
    ASSERT_FILE_EXISTS("lj.dump");
    ASSERT_FILE_EXISTS("rdf.first");

    // recomputes the RDF from the dump with 100 bins out to 5 sigma
    run_input("in.rdf.rerun");
    ASSERT_FILE_EXISTS("rdf.rerun");

    auto first = last_vector_block("rdf.first");
    auto rerun = last_vector_block("rdf.rerun");
    ASSERT_EQ(first.size(), 50);
    ASSERT_EQ(rerun.size(), 100);

    // the first 50 bins cover the same distances with the same 0.05 sigma
    // bin width, so the g(r) values must match apart from pairs that the
    // "%g" precision of the dump file moves across a bin boundary
    for (std::size_t i = 0; i < first.size(); ++i) {
        ASSERT_EQ(first[i].size(), 4);
        ASSERT_EQ(rerun[i].size(), 4);
        EXPECT_NEAR(first[i][1], rerun[i][1], 1.0e-12);
        EXPECT_NEAR(first[i][2], rerun[i][2], 0.05);
    }

    delete_file("lj.dump");
    delete_file("rdf.first");
    delete_file("rdf.rerun");
}

} // namespace LAMMPS_NS

EXAMPLE_TEST_MAIN()
