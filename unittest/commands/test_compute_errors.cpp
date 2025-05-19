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

#include "../testing/core.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class ComputeErrors : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ComputeErrors";
        LAMMPSTest::SetUp();
        command("atom_style atomic");
    }
};

TEST_F(ComputeErrors, PropertyAtom)
{

    std::vector<std::string> properties_not_available = {
        "mol", "q", "mux", "muy", "muz", "mu",
        "spx", "spy", "spz", "sp", "fmx", "fmy", "fmz",
        "nbonds",
        "radius", "diameter", "omegax", "omegay", "omegaz",
        "temperature", "heatflow",
        "angmomx", "angmomy", "angmomz",
        "quatw", "quati", "quatj", "quatk", "tqx", "tqy", "tqz"
    };

    for (const std::string& property : properties_not_available)
        TEST_FAILURE(".*ERROR: Compute property/atom " + property + " is not available.*",
            command("compute 1 all property/atom " + property););

    std::vector<std::string> properties_requires = {
        "shapex", "shapey", "shapez",
        "end1x", "end1y", "end1z", "end2x", "end2y", "end2z",
        "corner1x", "corner1y", "corner1z",
        "corner2x", "corner2y", "corner2z",
        "corner3x", "corner3y", "corner3z"
    };

    for (const std::string& property : properties_requires)
        TEST_FAILURE(".*ERROR: Compute property/atom " + property + " requires atom style.*",
            command("compute 1 all property/atom " + property););

    std::vector<std::string> properties_custom = {
        "i_name", "d_name", "i2_name[0]", "d2_name[0]",
    };

    for (const std::string& property : properties_custom)
        TEST_FAILURE(".*ERROR: Compute property/atom .* does not exist.*",
            command("compute 1 all property/atom " + property););

    std::vector<std::string> properties_others = {
        "vfrac", "s0", "espin", "eradius", "ervel", "erforce",
        "rho", "drho", "e", "de", "cv", "buckling"
    };

    for (const std::string& property : properties_others)
        TEST_FAILURE(".*ERROR: Invalid keyword " + property + " for atom style .* in compute property/atom command.*",
            command("compute 1 all property/atom " + property););

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
