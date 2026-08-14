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

#include "lammps.h"

#include "atom.h"
#include "force.h"
#include "info.h"
#include "pair.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

// write a two-species local-density potential file with the "noforce"
// potential so that the local densities themselves are isolated.
static void write_noforce_ldd(const char *fname)
{
    FILE *fp = fopen(fname, "w");
    fputs("A A indicator lucy   0.0 0.65 self yes potential noforce\n"
          "A B indicator dpd    0.0 0.65 self no  potential noforce\n"
          "B A indicator sphere 0.0 0.55 self no  potential noforce\n"
          "B B indicator lucy   0.0 0.65 self no  potential noforce\n",
          fp);
    fclose(fp);
}

class PairLDDTest : public LAMMPSTest {
protected:
    Atom *atom;
    Force *force;
    void SetUp() override
    {
        testbinary = "PairLDDTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite"};
        LAMMPSTest::SetUp();
        atom  = lmp->atom;
        force = lmp->force;
    }
    void TearDown() override { LAMMPSTest::TearDown(); }
};

// the local densities computed by pair ldd must be reachable through
// Pair::extract_peratom() (the path used by fix pair) and must match the
// hand-verified reference values for each indicator type.
TEST_F(PairLDDTest, local_density_values)
{
    if (!Info::has_package("BOCS")) GTEST_SKIP();
    write_noforce_ldd("test_ldd_noforce.ldd");

    BEGIN_HIDE_OUTPUT();
    command("units real");
    command("atom_style atomic");
    command("atom_modify map array");
    command("region my_box block 0 9 0 9 0 9");
    command("create_box 2 my_box");
    command("mass * 59.0448");
    command("pair_style ldd");
    command("create_atoms 1 single 1    0.1    0.75");
    command("create_atoms 2 single 1    0.1    0.65");
    command("create_atoms 1 single 1    0.1    0.1");
    command("create_atoms 2 single 1    0.54   0.65");
    command("create_atoms 2 single 1    8.9    0.65");
    command("pair_coeff * * test_ldd_noforce.ldd A B");
    command("neighbor 4.0 bin");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    int ncol  = -1;
    auto **ld = (double **) force->pair->extract_peratom("local_density", ncol);
    ASSERT_NE(ld, nullptr);
    ASSERT_EQ(ncol, 2);    // two species -> two columns (A = col 0, B = col 1)

    // column 0 = density of species A (type 1), column 1 = species B (type 2)
    const double ref[5][2] = {{7.606403744, 10.77780578}, {8.383039619, 5.631116274},
                              {7.606403744, 0.29205516},   {0.522157906, 0.777519395},
                              {4.864480249, 4.853815902}};
    for (int i = 0; i < 5; ++i) {
        EXPECT_NEAR(ld[i][0], ref[i][0], 1e-7);
        EXPECT_NEAR(ld[i][1], ref[i][1], 1e-7);
    }
    remove("test_ldd_noforce.ldd");
}

// multiple atom types may map to the same local-density species.  Splitting a
// species across two atom types (here A across types 1 and 3) must leave the
// computed densities unchanged, since the species assignment is identical.
TEST_F(PairLDDTest, same_species_invariance)
{
    if (!Info::has_package("BOCS")) GTEST_SKIP();
    write_noforce_ldd("test_ldd_noforce.ldd");

    // baseline: two types, mapped A B
    BEGIN_HIDE_OUTPUT();
    command("units real");
    command("atom_style atomic");
    command("atom_modify map array");
    command("region my_box block 0 9 0 9 0 9");
    command("create_box 2 my_box");
    command("mass * 59.0448");
    command("pair_style ldd");
    command("create_atoms 1 single 1    0.1    0.75");
    command("create_atoms 2 single 1    0.1    0.65");
    command("create_atoms 1 single 1    0.1    0.1");
    command("create_atoms 2 single 1    0.54   0.65");
    command("create_atoms 2 single 1    8.9    0.65");
    command("pair_coeff * * test_ldd_noforce.ldd A B");
    command("neighbor 4.0 bin");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    int ncol      = -1;
    auto **ld     = (double **) lmp->force->pair->extract_peratom("local_density", ncol);
    double a0_base = ld[0][0], b0_base = ld[0][1];    // species A/B density of atom 0

    // same geometry, but the type-1 atoms are relabeled as a third type that
    // also maps to species A; the mapping A B A must reproduce the densities
    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("atom_style atomic");
    command("atom_modify map array");
    command("region my_box block 0 9 0 9 0 9");
    command("create_box 3 my_box");
    command("mass * 59.0448");
    command("pair_style ldd");
    command("create_atoms 3 single 1    0.1    0.75");    // was type 1 -> species A
    command("create_atoms 2 single 1    0.1    0.65");
    command("create_atoms 3 single 1    0.1    0.1");      // was type 1 -> species A
    command("create_atoms 2 single 1    0.54   0.65");
    command("create_atoms 2 single 1    8.9    0.65");
    command("pair_coeff * * test_ldd_noforce.ldd A B A");
    command("neighbor 4.0 bin");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    ncol           = -1;
    auto **ld2     = (double **) lmp->force->pair->extract_peratom("local_density", ncol);
    double a0_split = ld2[0][0], b0_split = ld2[0][1];

    EXPECT_NEAR(a0_split, a0_base, 1e-10);
    EXPECT_NEAR(b0_split, b0_base, 1e-10);
    remove("test_ldd_noforce.ldd");
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
            if (arg == "-v") verbose = true;
        }
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
