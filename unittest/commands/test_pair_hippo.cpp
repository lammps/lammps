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
#include "comm.h"
#include "info.h"
#include "input.h"
#include "library.h"
#include "utils.h"
#include "variable.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// HIPPO water hexamer (the examples/amoeba system) without cutoff effects.
// Reference energies in kcal/mol from Tinker 8 (version 26.2) for the same
// parameters, which use the Pauli repulsion convention of Tinker 8.9 and later.

static constexpr double tinker_repulsion = 94.7840;
static constexpr double tinker_pair      = 94.7840 - 16.2997 - 85.5277 - 18.0250 - 22.2938;

static const char *key = "cutoff 10.0\n"
                         "taper 8.0\n"
                         "polar-eps 1.0e-10\n";

static std::string read_file(const std::string &path)
{
    std::ifstream in(path);
    std::stringstream buf;
    buf << in.rdbuf();
    return buf.str();
}

static void write_file(const std::string &path, const std::string &text)
{
    std::ofstream out(path);
    out << text;
}

class PairHippoTest : public LAMMPSTest {
protected:
    std::string prm;

    void SetUp() override
    {
        testbinary = "PairHippoTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite"};
        LAMMPSTest::SetUp();
        prm = read_file(STRINGIFY(TEST_INPUT_FOLDER) "/hippo_water.prm");
        write_file("test_hippo_water.key", key);
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        platform::unlink("test_hippo_water.key");
        platform::unlink("test_hippo_water.prm");
    }

    // set up the hexamer with the PRM file text and the pair style
    void setup_system(const std::string &prmtext, const std::string &style = "hippo")
    {
        write_file("test_hippo_water.prm", prmtext);
        command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
        command("units real");
        command("boundary s s s");
        command("atom_style amoeba");
        command("atom_modify map array sort 0 0.0");
        command("bond_style zero");
        command("angle_style zero");
        command("fix amtype all property/atom i_amtype ghost yes");
        command("fix extra all property/atom i_amgroup d_redID d_pval ghost yes");
        command("fix extra2 all property/atom i_polaxe d2_xyzaxis 3");
        command("read_data ${input_dir}/data.hippo_water fix amtype NULL \"Tinker Types\"");
        command("bond_coeff *");
        command("angle_coeff *");
        command("pair_style " + style);
        command("pair_coeff * * test_hippo_water.prm test_hippo_water.key");
        command("special_bonds lj/coul 0.5 0.5 0.5 one/five yes");
        command("compute pa all pair " + style);
        command("thermo_style custom step pe c_pa[2]");
    }

    double energy()
    {
        BEGIN_HIDE_OUTPUT();
        command("run 0 post no");
        END_HIDE_OUTPUT();
        return lammps_get_thermo(lmp, "pe");
    }

    double repulsion()
    {
        energy();
        return lmp->input->variable->compute_equal("c_pa[2]");
    }
};

// with the current Tinker convention the repulsion and the total nonbonded
// energy agree with Tinker

TEST_F(PairHippoTest, energy_matches_tinker)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    BEGIN_HIDE_OUTPUT();
    setup_system(prm);
    END_HIDE_OUTPUT();
    EXPECT_NEAR(repulsion(), tinker_repulsion, 1.0e-3);
    EXPECT_NEAR(energy(), tinker_pair, 1.0e-3);
}

// analytic forces are the gradient of the energy

TEST_F(PairHippoTest, forces_finite_differences)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    BEGIN_HIDE_OUTPUT();
    setup_system(prm);
    END_HIDE_OUTPUT();

    energy();
    const int natoms = lmp->atom->natoms;
    std::vector<double> fref(3 * natoms);
    for (int id = 1; id <= natoms; ++id) {
        const int i = lmp->atom->map(id);
        for (int k = 0; k < 3; ++k)
            fref[3 * (id - 1) + k] = lmp->atom->f[i][k];
    }

    const double delta = 1.0e-4;
    double maxerr      = 0.0;
    for (int id = 1; id <= natoms; ++id) {
        for (int k = 0; k < 3; ++k) {
            lmp->atom->x[lmp->atom->map(id)][k] += delta;
            const double eplus = energy();
            lmp->atom->x[lmp->atom->map(id)][k] -= 2.0 * delta;
            const double eminus = energy();
            lmp->atom->x[lmp->atom->map(id)][k] += delta;
            const double fd = -(eplus - eminus) / (2.0 * delta);
            maxerr          = std::max(maxerr, std::fabs(fd - fref[3 * (id - 1) + k]));
        }
    }
    EXPECT_LT(maxerr, 1.0e-3);
}

// repulsion sizes of Tinker before version 8.9 are 8 times smaller: such a
// file gives a warning, and a repulsion 64 times smaller

TEST_F(PairHippoTest, old_convention_warning)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    std::string old = prm;

    const std::vector<std::pair<std::string, std::string>> sizes = {{"23.0064", "2.8758"},
                                                                    {"16.4880", "2.0610"}};
    for (const auto &s : sizes) {
        const auto pos = old.find(s.first);
        ASSERT_NE(pos, std::string::npos);
        old.replace(pos, s.first.size(), s.second);
    }
    BEGIN_CAPTURE_OUTPUT();
    setup_system(old);
    auto output = END_CAPTURE_OUTPUT();
    EXPECT_THAT(output, testing::HasSubstr("convention of Tinker before version 8.9"));
    EXPECT_NEAR(repulsion(), tinker_repulsion / 64.0, 1.0e-3);
}

// charge flux is not supported and must not be ignored

TEST_F(PairHippoTest, charge_flux_error)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    std::string cflux = prm + "\n\n      ##########################\n"
                              "      ##                      ##\n"
                              "      ##  Charge Flux Parameters  ##\n"
                              "      ##                      ##\n"
                              "      ##########################\n\n\n"
                              "bndcflux      1    2           -0.027265\n";
    TEST_FAILURE(".*ERROR: HIPPO PRM file line .*: charge flux \\(bndcflux\\) is not supported.*",
                 setup_system(cflux););
}

// d-equals-p is accepted for HIPPO and rejected for AMOEBA

TEST_F(PairHippoTest, d_equals_p)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    std::string text = prm;

    const std::string anchor = "polarization            MUTUAL\n";
    const auto pos           = text.find(anchor);
    ASSERT_NE(pos, std::string::npos);
    text.insert(pos + anchor.size(), "d-equals-p\n");
    BEGIN_HIDE_OUTPUT();
    setup_system(text);
    END_HIDE_OUTPUT();
    EXPECT_NEAR(repulsion(), tinker_repulsion, 1.0e-3);

    BEGIN_HIDE_OUTPUT();
    command("clear");
    END_HIDE_OUTPUT();
    std::string amoeba = read_file(STRINGIFY(TEST_INPUT_FOLDER) "/amoeba_water.prm");
    const auto apos    = amoeba.find(anchor);
    ASSERT_NE(apos, std::string::npos);
    amoeba.insert(apos + anchor.size(), "d-equals-p\n");
    TEST_FAILURE(".*ERROR: AMOEBA PRM file keyword d-equals-p on line .* is not supported.*",
                 setup_system(amoeba, "amoeba"););
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
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
