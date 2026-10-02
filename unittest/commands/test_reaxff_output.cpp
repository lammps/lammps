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
#include "../testing/utils.h"

#include "atom.h"
#include "info.h"
#include "lammps.h"
#include "library.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

// the fixes reaxff/bonds and reaxff/species of the REAXFF package do not
// contribute to forces or energies but write the bond topology and the
// molecular species detected by the reactive force field to a file.  the
// force-style YAML test drivers therefore cannot cover them.  the tests
// below build a small gas phase system of four water molecules and two
// hydrogen molecules, whose bond topology and species composition are known
// in advance, run a few steps, and parse the two output files back.  the
// ReaxFFOutputKokkos ctest entry runs the very same test bodies with
// LAMMPS_ACCELERATOR_ARGS="-k on t 1 -sf kk" so that the KOKKOS versions of
// the two fixes are exercised as well.

class ReaxFFOutputTest : public LAMMPSTest {
protected:
    // geometry of the molecules placed into the box.  the oxygen of each water
    // molecule is created first, so the atom-IDs are 1-3, 4-6, 7-9, and 10-12
    // for the four water molecules and 13-14 and 15-16 for the two hydrogen
    // molecules.  the molecules are 8 angstrom or more apart, which is well
    // beyond any bond order cutoff, so they stay individual molecules.
    static constexpr int num_water    = 4;
    static constexpr int num_hydrogen = 2;

    void SetUp() override
    {
        testbinary = "ReaxFFOutputTest";
        LAMMPSTest::SetUp();
        if (!info->has_style("pair", "reaxff")) return;

        const double water[3][3] = {{0.0, 0.0, 0.0}, {0.9572, 0.0, 0.0}, {-0.2400, 0.9266, 0.0}};
        const double hydro[2][3] = {{0.0, 0.0, 0.0}, {0.7414, 0.0, 0.0}};
        const double wpos[num_water][3] = {
            {3.0, 3.0, 3.0}, {9.0, 3.0, 9.0}, {3.0, 9.0, 9.0}, {9.0, 9.0, 3.0}};
        const double hpos[num_hydrogen][3] = {{3.0, 3.0, 9.0}, {9.0, 9.0, 9.0}};

        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("atom_style charge");
        command("atom_modify map array");
        command("boundary p p p");
        command("region box block 0 12 0 12 0 12");
        command("create_box 2 box");
        command("mass 1 1.008");
        command("mass 2 15.999");

        // atom type 1 is hydrogen, atom type 2 is oxygen
        for (int m = 0; m < num_water; ++m)
            for (int a = 0; a < 3; ++a)
                command(fmt::format("create_atoms {} single {:.6f} {:.6f} {:.6f} units box",
                                    (a == 0) ? 2 : 1, wpos[m][0] + water[a][0],
                                    wpos[m][1] + water[a][1], wpos[m][2] + water[a][2]));
        for (int m = 0; m < num_hydrogen; ++m)
            for (int a = 0; a < 2; ++a)
                command(fmt::format("create_atoms 1 single {:.6f} {:.6f} {:.6f} units box",
                                    hpos[m][0] + hydro[a][0], hpos[m][1] + hydro[a][1],
                                    hpos[m][2] + hydro[a][2]));

        command("set type * charge 0.0");
        command("pair_style reaxff NULL");
        command("pair_coeff * * ffield.reax.mattsson H O");
        command("fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff");
        command("neighbor 1.0 bin");
        command("neigh_modify every 1 delay 0 check no");
        command("timestep 0.1");
        command("fix nve all nve");
        END_HIDE_OUTPUT();
    }

    // number of atoms in a chemical formula like "H2O" or "C2H4".
    // returns -1 if the string is not a valid formula.
    static int formula_atoms(const std::string &formula)
    {
        int total       = 0;
        std::size_t idx = 0;
        while (idx < formula.size()) {
            if ((formula[idx] < 'A') || (formula[idx] > 'Z')) return -1;
            ++idx;
            while ((idx < formula.size()) && (formula[idx] >= 'a') && (formula[idx] <= 'z'))
                ++idx;
            int count = 0;
            while ((idx < formula.size()) && (formula[idx] >= '0') && (formula[idx] <= '9')) {
                count = count * 10 + (formula[idx] - '0');
                ++idx;
            }
            if (count == 0) count = 1;
            total += count;
        }
        return total;
    }

    // one atom record of the fix reaxff/bonds output file:
    // id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q
    struct BondRecord {
        int type;
        std::vector<int> neighid;
        std::vector<double> bondorder;
        double abo;
        double nlp;
        double q;
    };

    struct BondBlock {
        bigint timestep;
        bigint natoms;
        std::map<int, BondRecord> records;
    };
};

TEST_F(ReaxFFOutputTest, Bonds)
{
    if (!info->has_style("pair", "reaxff")) GTEST_SKIP();
    if (!info->has_style("fix", "reaxff/bonds")) GTEST_SKIP();

    const std::string outfile = "test_reaxff_bonds.tmp";
    delete_file(outfile);

    BEGIN_HIDE_OUTPUT();
    command("fix bonds all reaxff/bonds 5 " + outfile);
    command("run 10 post no");
    END_HIDE_OUTPUT();

    const int natoms = (int)lmp->atom->natoms;
    ASSERT_EQ(natoms, 3 * num_water + 2 * num_hydrogen);

    ASSERT_FILE_EXISTS(outfile);
    auto lines = read_lines(outfile);
    ASSERT_GT((int)lines.size(), 0);

    // split the file into the individual snapshots.  everything that starts
    // with a '#' is a comment, only "# Timestep" and "# Number of particles"
    // carry information that is checked here.
    std::vector<BondBlock> blocks;
    for (const auto &line : lines) {
        auto words = utils::split_words(line);
        if (words.empty()) continue;
        if (words[0] == "#") {
            if ((words.size() > 2) && (words[1] == "Timestep")) {
                BondBlock block;
                block.timestep = std::stoll(words[2]);
                block.natoms   = 0;
                blocks.push_back(block);
            } else if ((words.size() > 4) && (words[1] == "Number") && (words[2] == "of")) {
                ASSERT_FALSE(blocks.empty());
                blocks.back().natoms = std::stoll(words[4]);
            }
            continue;
        }

        ASSERT_FALSE(blocks.empty());
        // id type nb, nb neighbor IDs, molecule ID, nb bond orders, abo, nlp, q
        ASSERT_GE((int)words.size(), 3);
        const int id = std::stoi(words[0]);
        BondRecord rec;
        rec.type     = std::stoi(words[1]);
        const int nb = std::stoi(words[2]);
        ASSERT_EQ((int)words.size(), 7 + 2 * nb);
        for (int k = 0; k < nb; ++k)
            rec.neighid.push_back(std::stoi(words[3 + k]));
        for (int k = 0; k < nb; ++k)
            rec.bondorder.push_back(std::stod(words[4 + nb + k]));
        rec.abo = std::stod(words[4 + 2 * nb]);
        rec.nlp = std::stod(words[5 + 2 * nb]);
        rec.q   = std::stod(words[6 + 2 * nb]);
        ASSERT_EQ((int)blocks.back().records.count(id), 0);
        blocks.back().records[id] = rec;
    }

    // fix reaxff/bonds with nevery 5 writes at the setup of the run and then
    // every 5 steps, i.e. for timesteps 0, 5, and 10
    ASSERT_EQ((int)blocks.size(), 3);
    EXPECT_EQ(blocks[0].timestep, 0);
    EXPECT_EQ(blocks[1].timestep, 5);
    EXPECT_EQ(blocks[2].timestep, 10);

    for (const auto &block : blocks) {
        // the header must report the number of atoms and there must be exactly
        // one record for each atom, with the IDs covering the full range
        EXPECT_EQ(block.natoms, natoms);
        ASSERT_EQ((int)block.records.size(), natoms);
        for (int id = 1; id <= natoms; ++id)
            EXPECT_EQ((int)block.records.count(id), 1);

        int nbonds   = 0;
        double sumq  = 0.0;
        double maxbo = 0.0;
        for (const auto &item : block.records) {
            const int id        = item.first;
            const BondRecord &r = item.second;
            const int nb        = (int)r.bondorder.size();

            // hydrogen (type 1) has exactly one bond, oxygen (type 2) has two
            EXPECT_EQ(nb, (r.type == 2) ? 2 : 1);
            EXPECT_EQ((int)r.neighid.size(), nb);

            double sumbo = 0.0;
            for (int k = 0; k < nb; ++k) {
                // bond orders must be physically plausible: larger than the
                // 0.3 bond graph cutoff that selects them and well below the
                // maximum of 4 that a single covalent bond can reach
                EXPECT_GT(r.bondorder[k], 0.3);
                EXPECT_LT(r.bondorder[k], 4.0);
                sumbo += r.bondorder[k];
                maxbo = std::max(maxbo, r.bondorder[k]);

                // the bond must be listed by the bonded partner as well
                const int jd = r.neighid[k];
                EXPECT_GE(jd, 1);
                EXPECT_LE(jd, natoms);
                ASSERT_EQ((int)block.records.count(jd), 1);
                const auto &partner = block.records.at(jd);
                EXPECT_NE(std::find(partner.neighid.begin(), partner.neighid.end(), id),
                          partner.neighid.end());
            }
            nbonds += nb;

            // the sum of the bond orders of the listed bonds must not exceed
            // the total bond order of the atom (up to the rounding of the
            // printed values) and the difference, i.e. the bonds below the
            // cutoff, must be small for these isolated molecules
            EXPECT_GT(r.abo - sumbo, -0.005);
            EXPECT_LT(r.abo - sumbo, 0.1);

            // oxygen has two lone pairs, hydrogen has none
            if (r.type == 2) {
                EXPECT_NEAR(r.nlp, 2.0, 0.1);
            } else {
                EXPECT_NEAR(r.nlp, 0.0, 0.1);
            }
            sumq += r.q;
        }

        // every hydrogen contributes one and every oxygen two bonds, and each
        // bond is listed by both of its atoms
        EXPECT_EQ(nbonds % 2, 0);
        EXPECT_EQ(nbonds / 2, 2 * num_water + num_hydrogen);

        // the bond orders must not be all zero
        EXPECT_GT(maxbo, 0.5);

        // charge equilibration conserves the total charge of the neutral system
        EXPECT_NEAR(sumq, 0.0, 0.02);
    }

    // representative values of the undistorted initial geometry: the O-H bond
    // of a water molecule, the total bond order and charge of its oxygen, and
    // the stronger H-H bond of a hydrogen molecule
    const auto &first = blocks[0].records;
    EXPECT_NEAR(first.at(1).bondorder[0], 0.936, 0.01);
    EXPECT_NEAR(first.at(1).abo, 1.873, 0.01);
    EXPECT_NEAR(first.at(1).q, -0.653, 0.01);
    EXPECT_NEAR(first.at(2).bondorder[0], 0.936, 0.01);
    EXPECT_NEAR(first.at(13).bondorder[0], 0.966, 0.01);
    EXPECT_GT(first.at(13).bondorder[0], first.at(1).bondorder[0]);

    delete_file(outfile);
}

TEST_F(ReaxFFOutputTest, Species)
{
    if (!info->has_style("pair", "reaxff")) GTEST_SKIP();
    if (!info->has_style("fix", "reaxff/species")) GTEST_SKIP();

    const std::string outfile = "test_reaxff_species.tmp";
    delete_file(outfile);

    BEGIN_HIDE_OUTPUT();
    command("fix species all reaxff/species 1 1 5 " + outfile);
    command("run 10 post no");
    END_HIDE_OUTPUT();

    const int natoms = (int)lmp->atom->natoms;
    ASSERT_EQ(natoms, 3 * num_water + 2 * num_hydrogen);

    ASSERT_FILE_EXISTS(outfile);
    auto lines = read_lines(outfile);

    // each snapshot consists of a header line naming the species found and a
    // data line with the timestep, the number of molecules, the number of
    // species, and the number of molecules of each species
    ASSERT_EQ((int)lines.size(), 6);

    const bigint expected_steps[3] = {0, 5, 10};
    for (int n = 0; n < 3; ++n) {
        auto head = utils::split_words(lines[2 * n]);
        auto data = utils::split_words(lines[2 * n + 1]);

        ASSERT_GE((int)head.size(), 4);
        EXPECT_EQ(head[0], "#");
        EXPECT_EQ(head[1], "Timestep");
        EXPECT_EQ(head[2], "No_Moles");
        EXPECT_EQ(head[3], "No_Specs");

        ASSERT_GE((int)data.size(), 3);
        EXPECT_EQ(std::stoll(data[0]), expected_steps[n]);
        const int nmole = std::stoi(data[1]);
        const int nspec = std::stoi(data[2]);

        // one column with the number of molecules for each species
        ASSERT_EQ((int)head.size(), 4 + nspec);
        ASSERT_EQ((int)data.size(), 3 + nspec);

        std::map<std::string, int> species;
        int summole = 0;
        int sumatom = 0;
        for (int k = 0; k < nspec; ++k) {
            const std::string name = head[4 + k];
            const int count        = std::stoi(data[3 + k]);
            EXPECT_GT(count, 0);
            const int size = formula_atoms(name);
            EXPECT_GT(size, 0) << "not a valid chemical formula: " << name;
            summole += count;
            sumatom += count * size;
            species[name] = count;
        }

        // the molecule counts of the individual species add up to the total
        // number of molecules and the atoms of all molecules to all atoms
        EXPECT_EQ(summole, nmole);
        EXPECT_EQ(sumatom, natoms);

        // the four water molecules and the two hydrogen molecules must be
        // recognized as such and stay intact for the duration of the run
        EXPECT_EQ(nspec, 2);
        EXPECT_EQ(nmole, num_water + num_hydrogen);
        EXPECT_EQ(species["H2O"], num_water);
        EXPECT_EQ(species["H2"], num_hydrogen);
    }

    delete_file(outfile);
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

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
