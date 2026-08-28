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
// there prescribes and check the claims it makes.  The regression tests skip
// the three dependent inputs (in.rerun, in.read_dump, in.rdf.rerun): they
// read the lj.dump file that only a full run of in.first or in.rdf.first
// writes, and since both writers use the same file name, running the folder
// with independent workers makes the results depend on scheduling.

#include "../testing/core.h"
#include "../testing/utils.h"

#include "output.h"
#include "thermo.h"
#include "update.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class RerunExamplesTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "RerunExamplesTest";
        LAMMPSTest::SetUp();
    }

    void run_input(const std::string &script)
    {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        command("include \"" STRINGIFY(TEST_EXAMPLES_FOLDER) "/" + script + "\"");
        END_HIDE_OUTPUT();
    }

    double thermo_value(const std::string &keyword)
    {
        double value = 0.0;
        lmp->output->thermo->evaluate_keyword(keyword, &value);
        return value;
    }
};

// last output block of a fix ave/time "mode vector" file, as rows of columns
static std::vector<std::vector<double>> last_vector_block(const std::string &filename)
{
    std::vector<std::vector<double>> block;
    std::ifstream data(filename);
    if (!data.is_open()) return block;

    std::string line;
    std::size_t nrows = 0;
    while (std::getline(data, line)) {
        auto words = utils::split_words(line);
        if (words.empty() || (words[0][0] == '#')) continue;
        if (words.size() == 2) { // "timestep number-of-rows" starts a block
            nrows = std::stoul(words[1]);
            block.clear();
        } else {
            std::vector<double> row;
            for (const auto &word : words)
                row.push_back(std::stod(word));
            block.push_back(row);
        }
    }
    EXPECT_EQ(block.size(), nrows);
    return block;
}

TEST_F(RerunExamplesTest, first_rerun_read_dump)
{
    // writes the snapshots 0 to 1000 of a 32000 atom LJ melt to lj.dump
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

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
