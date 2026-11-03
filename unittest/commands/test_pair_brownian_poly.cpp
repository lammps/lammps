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

// the pairwise random forces of pair style brownian/poly must be equal and opposite
// and must not depend on the number of MPI ranks or OpenMP threads

#include "lammps.h"

#include "atom.h"
#include "info.h"
#include "input.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
class BrownianPolyTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "BrownianPolyTest";
        LAMMPSTest::SetUp();
        if (!info->has_style("pair", "brownian/poly")) GTEST_SKIP() << "COLLOID package not enabled";
    }

    void make_box()
    {
        BEGIN_HIDE_OUTPUT();
        command("units lj");
        command("atom_style sphere");
        command("atom_modify map array");
        command("boundary p p p");
        command("newton off");
        command("comm_modify mode single vel yes");
        command("region box block -20 20 -20 20 -20 20");
        command("create_box 1 box");
        END_HIDE_OUTPUT();
    }

    // two spheres of unequal size close to contact across the periodic boundary
    void make_pair()
    {
        BEGIN_HIDE_OUTPUT();
        command("create_atoms 1 single 19.5 0.3 -0.4");
        command("create_atoms 1 single -17.0 2.1 1.9");
        command("set atom 1 diameter 2.0 density 1.0");
        command("set atom 2 diameter 6.0 density 1.0");
        END_HIDE_OUTPUT();
    }

    // per-atom vector of the atom with ID "id" (on any MPI rank)
    void get_atom_vec(double **vec, tagint id, double *val)
    {
        double tmp[3] = {0.0, 0.0, 0.0};
        const int i   = lmp->atom->map(id);
        if ((i >= 0) && (i < lmp->atom->nlocal))
            for (int k = 0; k < 3; ++k) tmp[k] = vec[i][k];
        MPI_Allreduce(tmp, val, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    // net force on all atoms
    void net_force(double *fnet)
    {
        double tmp[3] = {0.0, 0.0, 0.0};
        double **f    = lmp->atom->f;
        for (int i = 0; i < lmp->atom->nlocal; ++i)
            for (int k = 0; k < 3; ++k) tmp[k] += f[i][k];
        MPI_Allreduce(tmp, fnet, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
};

// forces must be equal and opposite at every timestep
TEST_F(BrownianPolyTest, Newton3rdLawTwoSpheres)
{
    make_box();
    make_pair();
    BEGIN_HIDE_OUTPUT();
    // mu flaglog flagfld cutinner cutoff t_target seed flagHI flagVF
    command("pair_style brownian/poly 1.0 1 0 4.1 6.0 1.0 5827 1 0");
    command("pair_coeff * *");
    END_HIDE_OUTPUT();

    double f1[3], f2[3];
    for (int step = 0; step < 5; ++step) {
        BEGIN_HIDE_OUTPUT();
        command("run 1 post no");
        END_HIDE_OUTPUT();
        get_atom_vec(lmp->atom->f, 1, f1);
        get_atom_vec(lmp->atom->f, 2, f2);
        EXPECT_GT(fabs(f1[0]), 1.0e-3);
        for (int k = 0; k < 3; ++k)
            EXPECT_NEAR(f1[k] + f2[k], 0.0, 1.0e-12) << "step " << step << " component " << k;
    }
}

// same with the centers inside cutinner, so that the minimum gap is used
TEST_F(BrownianPolyTest, Newton3rdLawClampedGap)
{
    make_box();
    make_pair();
    BEGIN_HIDE_OUTPUT();
    command("pair_style brownian/poly 1.0 1 0 5.5 6.0 1.0 91237 1 0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    double f1[3], f2[3];
    get_atom_vec(lmp->atom->f, 1, f1);
    get_atom_vec(lmp->atom->f, 2, f2);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(f1[k] + f2[k], 0.0, 1.0e-12) << "component " << k;
}

// cluster of polydisperse spheres near contact: the net force must vanish.  the
// cutoff excludes the diagonal pairs, where the lubrication resistances are negative
TEST_F(BrownianPolyTest, Newton3rdLawCluster)
{
    make_box();
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("create_atoms 1 single 2.4 0.0 0.0");
    command("create_atoms 1 single 0.0 2.4 0.0");
    command("create_atoms 1 single 2.4 2.4 0.0");
    command("set atom 1 diameter 2.0 density 1.0");
    command("set atom 2 diameter 2.2 density 1.0");
    command("set atom 3 diameter 2.1 density 1.0");
    command("set atom 4 diameter 2.3 density 1.0");
    command("pair_style brownian/poly 1.0 1 0 0.5 3.0 1.0 40591 1 0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    double fnet[3];
    net_force(fnet);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(fnet[k], 0.0, 1.0e-12) << "component " << k;
}

// the random numbers depend only on the atom IDs, the timestep, and the seed.
// compare to reference data, so that a dependence on the number of MPI ranks
// or OpenMP threads is detected when the test is run in parallel
TEST_F(BrownianPolyTest, ReferenceForces)
{
    make_box();
    make_pair();
    BEGIN_HIDE_OUTPUT();
    command("pair_style brownian/poly 1.0 1 0 4.1 6.0 1.0 5827 1 0");
    command("pair_coeff * *");
    command("reset_timestep 100");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const double f1_ref[3] = {-131.709970930346, -44.1376642992628, 53.3527990497969};
    const double t1_ref[3] = {-24.6292226024889, -15.9034629414138, 49.9254401753286};
    const double t2_ref[3] = {308.637940484353, -562.772959289939, -29.2354195536285};
    double f1[3], t1[3], t2[3];
    get_atom_vec(lmp->atom->f, 1, f1);
    get_atom_vec(lmp->atom->torque, 1, t1);
    get_atom_vec(lmp->atom->torque, 2, t2);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(f1[k], f1_ref[k], 1.0e-10 * fabs(f1_ref[k])) << "force component " << k;
        EXPECT_NEAR(t1[k], t1_ref[k], 1.0e-10 * fabs(t1_ref[k])) << "torque component " << k;
        EXPECT_NEAR(t2[k], t2_ref[k], 1.0e-10 * fabs(t2_ref[k])) << "torque component " << k;
    }
}

// the volume fraction correction must be the same whether it is computed
// once during setup or on every step because of fix deform
TEST_F(BrownianPolyTest, VolumeFractionDeform)
{
    double fref[3], fdef[3];
    for (int pass = 0; pass < 2; ++pass) {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        END_HIDE_OUTPUT();
        make_box();
        make_pair();
        BEGIN_HIDE_OUTPUT();
        // isotropic terms only, with volume fraction correction
        command("pair_style brownian/poly 1.0 1 1 4.1 6.0 1.0 5827 0 1");
        command("pair_coeff * *");
        if (pass == 1) command("fix 1 all deform 1 x final -20.0 20.0 remap v");
        command("run 0 post no");
        END_HIDE_OUTPUT();
        get_atom_vec(lmp->atom->f, 1, (pass == 0) ? fref : fdef);
    }
    EXPECT_GT(fabs(fref[0]), 1.0e-3);
    for (int k = 0; k < 3; ++k) EXPECT_DOUBLE_EQ(fdef[k], fref[k]) << "component " << k;
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
