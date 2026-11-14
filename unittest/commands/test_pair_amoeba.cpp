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
#include "domain.h"
#include "force.h"
#include "info.h"
#include "library.h"
#include "pair.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// AMOEBA water hexamer in a 10 Angstrom periodic box with PME.  With no
// ewald-alpha or pewald-alpha keyword, both Ewald coefficients are the same,
// so polar_kspace() re-uses the permanent multipole potential computed by
// multipole_kspace().  The second key file uses different coefficients.

static const char *key_equal = "ewald\n"
                               "ewald-cutoff 6.0\n"
                               "vdw-cutoff 6.0\n"
                               "vdw-taper 5.0\n"
                               "pme-grid 16 16 16\n"
                               "polar-eps 1.0e-10\n";

static const char *key_unequal = "ewald\n"
                                 "ewald-cutoff 6.0\n"
                                 "ewald-alpha 0.55\n"
                                 "pewald-alpha 0.65\n"
                                 "vdw-cutoff 6.0\n"
                                 "vdw-taper 5.0\n"
                                 "pme-grid 16 16 16\n"
                                 "polar-eps 1.0e-10\n";

// finite difference steps and tolerances. Single precision FFTs make the
// energy noisier, so use larger steps and tolerances there. Without the fix
// for equal Ewald coefficients, the errors are 1.2 kcal/mol/Angstrom in the
// forces and 1.9 kcal/mol in the virial.
#if defined(FFT_SINGLE)
static constexpr double force_delta = 1.0e-3;
static constexpr double force_tol   = 5.0e-3;
static constexpr double strain      = 1.0e-4;
static constexpr double virial_tol  = 2.0e-2;
static constexpr double predict_tol = 1.0e-3;
#else
static constexpr double force_delta = 1.0e-4;
static constexpr double force_tol   = 1.0e-3;
static constexpr double strain      = 1.0e-5;
static constexpr double virial_tol  = 1.0e-4;
static constexpr double predict_tol = 1.0e-4;
#endif

class PairAmoebaTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "PairAmoebaTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite"};
        LAMMPSTest::SetUp();
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        platform::unlink("test_amoeba_water.key");
    }

    void setup_system(const char *key)
    {
        FILE *fp = fopen("test_amoeba_water.key", "w");
        fputs(key, fp);
        fclose(fp);

        BEGIN_HIDE_OUTPUT();
        command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
        command("units real");
        command("boundary p p p");
        command("atom_style amoeba");
        command("atom_modify map array sort 0 0.0");
        command("bond_style zero");
        command("angle_style zero");
        command("fix amtype all property/atom i_amtype ghost yes");
        command("fix extra all property/atom i_amgroup d_redID d_pval ghost yes");
        command("fix extra2 all property/atom i_polaxe d2_xyzaxis 3");
        command("read_data ${input_dir}/data.amoeba_water fix amtype NULL \"Tinker Types\"");
        command("bond_coeff *");
        command("angle_coeff *");
        command("pair_style amoeba");
        command("pair_coeff * * ${input_dir}/amoeba_water.prm test_amoeba_water.key");
        command("special_bonds lj/coul 0.5 0.5 0.5 one/five yes");
        command("thermo_style custom step pe press");
        command("run 0 post no");
        END_HIDE_OUTPUT();
    }

    double energy()
    {
        BEGIN_HIDE_OUTPUT();
        command("run 0 post no");
        END_HIDE_OUTPUT();
        return lammps_get_thermo(lmp, "pe");
    }

    // largest difference between the analytic forces and central finite
    // differences of the potential energy, over all atoms and directions

    double max_force_error(double delta)
    {
        energy();
        const int natoms = lmp->atom->natoms;
        std::vector<double> fref(3 * natoms);
        for (int id = 1; id <= natoms; ++id) {
            const int i = lmp->atom->map(id);
            for (int k = 0; k < 3; ++k)
                fref[3 * (id - 1) + k] = lmp->atom->f[i][k];
        }

        double maxerr = 0.0;
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
        return maxerr;
    }

    // largest difference between the diagonal of the pair virial and the
    // derivative of the energy with respect to a uniform strain along x, y, z

    double max_virial_error(double strain)
    {
        energy();
        double vref[3];
        for (int k = 0; k < 3; ++k)
            vref[k] = lmp->force->pair->virial[k];

        const char *dim[3] = {"x", "y", "z"};
        double maxerr      = 0.0;
        for (int k = 0; k < 3; ++k) {
            BEGIN_HIDE_OUTPUT();
            command(fmt::format("change_box all {} scale {:.15g} remap units box", dim[k],
                                1.0 + strain));
            END_HIDE_OUTPUT();
            const double eplus = energy();
            BEGIN_HIDE_OUTPUT();
            command(fmt::format("change_box all {} scale {:.15g} remap units box", dim[k],
                                (1.0 - strain) / (1.0 + strain)));
            END_HIDE_OUTPUT();
            const double eminus = energy();
            BEGIN_HIDE_OUTPUT();
            command(fmt::format("change_box all {} scale {:.15g} remap units box", dim[k],
                                1.0 / (1.0 - strain)));
            END_HIDE_OUTPUT();
            const double dedstrain = (eplus - eminus) / (2.0 * strain);
            maxerr                 = std::max(maxerr, std::fabs(dedstrain + vref[k]));
        }
        return maxerr;
    }

    // energy after a single induced dipole iteration started from the
    // predicted dipoles, minus the converged energy. Without time integration
    // the atoms do not move, so the history of the predictor fills with the
    // converged dipoles of this geometry, and a correct prediction of the
    // next dipoles is those dipoles again.

    double predictor_error(const char *predictor)
    {
        const std::string key = std::string(key_unequal) + "polar-predict " + predictor + "\n";
        setup_system(key.c_str());
        const double converged = energy();
        BEGIN_HIDE_OUTPUT();
        command("run 20 post no");
        END_HIDE_OUTPUT();

        FILE *fp = fopen("test_amoeba_water.key", "w");
        fputs((key + "polar-iter 1\n").c_str(), fp);
        fclose(fp);
        BEGIN_HIDE_OUTPUT();
        command("pair_coeff * * ${input_dir}/amoeba_water.prm test_amoeba_water.key");
        END_HIDE_OUTPUT();
        return energy() - converged;
    }
};

// with equal Ewald coefficients, udirect1() overwrote the multipole potential
// that polar_kspace() re-uses, so the polarization forces were not the
// gradient of the energy (errors of up to 1.2 kcal/mol/Angstrom here).

TEST_F(PairAmoebaTest, forces_equal_alpha)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    setup_system(key_equal);
    EXPECT_LT(max_force_error(force_delta), force_tol);
}

// ... and the virial was wrong (by 1.9 kcal/mol here).

TEST_F(PairAmoebaTest, virial_equal_alpha)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    setup_system(key_equal);
    EXPECT_LT(max_virial_error(strain), virial_tol);
}

TEST_F(PairAmoebaTest, forces_unequal_alpha)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    setup_system(key_unequal);
    EXPECT_LT(max_force_error(force_delta), force_tol);
}

TEST_F(PairAmoebaTest, virial_unequal_alpha)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    setup_system(key_unequal);
    EXPECT_LT(max_virial_error(strain), virial_tol);
}

// the coefficients of the Gear and ASPC predictors were never set, so the
// predicted dipoles came from uninitialized memory. Here one iteration from
// them was off by 0.92 kcal/mol; one iteration from the direct field, without
// a predictor, is off by 0.04 kcal/mol.

TEST_F(PairAmoebaTest, predictor_gear)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    EXPECT_NEAR(predictor_error("gear"), 0.0, predict_tol);
}

TEST_F(PairAmoebaTest, predictor_aspc)
{
    if (!Info::has_package("AMOEBA")) GTEST_SKIP();
    if (lmp->comm->nprocs > 1) GTEST_SKIP();
    EXPECT_NEAR(predictor_error("aspc"), 0.0, predict_tol);
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
