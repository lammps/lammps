/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#include "lammps.h"
#include "atom.h"
#include "compute.h"
#include "library.h"
#include "fix.h"
#include "info.h"
#include "modify.h"

#include "platform.h"

#include "../testing/core.h"
#include "gtest/gtest.h"

#include <mpi.h>

bool verbose = false;

namespace LAMMPS_NS {

class FixUVTTest : public LAMMPSTest {
protected:
    void require_extra_fix()
    {
        if (!Info::has_package("EXTRA-FIX")) GTEST_SKIP();
    }

    void setup_quadratic_system()
    {
        require_extra_fix();
        command("units lj");
        command("atom_style atomic");
        command("boundary p p p");
        command("lattice fcc 0.8442");
        command("region box block 0 3 0 3 0 3");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("mass 1 1.0");
        command("pair_style zero 2.5");
        command("pair_coeff * *");
        command("neighbor 0.3 bin");
        command("neigh_modify every 1 delay 0 check yes");
        command("velocity all create 1.0 492845 mom yes rot no dist gaussian");
        command("timestep 0.005");
    }

    void setup_quadratic_fix()
    {
        command("variable k_quad equal 5.0");
        command("variable N0_quad equal 1.0");
        command("variable dEdN equal v_k_quad*(f_cp[13]-v_N0_quad)");
        command("fix cp all uvt temp 1.0 1.0 0.5 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
    }

    double fix_value(const char *id, int index)
    {
        void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
        double value = *(double *) ptr;
        lammps_free(ptr);
        return value;
    }
};

// A force-stage provider detects reads before setup and extra reads in the
// integrator.  Its callback order must precede FixUVT's post_force callback.
class FixDEDNProvider : public Fix {
public:
    int reads = 0;
    bool ready = false;
    double value = 0.0;
    FixDEDNProvider(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
    {
        scalar_flag = 1;
        global_freq = 1;
        extscalar = 0;
    }
    int setmask() override { return FixConst::POST_FORCE | FixConst::POST_FORCE_RESPA; }
    void setup(int vflag) override { post_force(vflag); }
    void post_force(int) override
    {
        auto *cp = modify->get_fix_by_id("cp");
        int dim = 0;
        value = 5.0 * (*static_cast<double *>(cp->extract("ne", dim)) - 1.0);
        ready = true;
    }
    void post_force_respa(int vflag, int, int) override { post_force(vflag); }
    double compute_scalar() override
    {
        EXPECT_TRUE(ready);
        ++reads;
        return value;
    }
};

TEST_F(FixUVTTest, DerivativeIsInitializedAndUpdatedAfterCoordinateDrift)
{
    setup_quadratic_system();
    setup_quadratic_fix();
    command("run 0 post no");
    EXPECT_NEAR(fix_value("cp", 14), 4.0, 1.0e-12);
    command("run 1 post no");
    const double ne = fix_value("cp", 12);
    EXPECT_LT(ne, 1.8);
    EXPECT_NEAR(fix_value("cp", 14), 5.0 * (ne - 1.0), 1.0e-12);
}

TEST_F(FixUVTTest, ForceStageProviderIsReadOncePerStep)
{
    setup_quadratic_system();
    (*lmp->modify->fix_map)["dedn/provider"] =
        [](LAMMPS *lmp, int narg, char **arg) -> Fix * {
            return new FixDEDNProvider(lmp, narg, arg);
        };
    command("fix provider all dedn/provider");
    command("fix cp all uvt temp 1 1 0.5 mu 2 2 0.5 ne 1.8 dedn f_provider");
    command("run 10 post no");
    auto *provider = static_cast<FixDEDNProvider *>(lmp->modify->get_fix_by_id("provider"));
    EXPECT_EQ(provider->reads, 11);
    EXPECT_NEAR(fix_value("cp", 14), 5.0 * (fix_value("cp", 12) - 1.0), 1.0e-12);
    command("run 10 post no");
    EXPECT_EQ(provider->reads, 22);
}

class FixEarlyCompute : public Fix {
public:
    FixEarlyCompute(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {}
    int setmask() override { return FixConst::INITIAL_INTEGRATE; }
    void initial_integrate(int) override
    {
        auto *compute = modify->get_compute_by_id("position");
        compute->compute_vector();
        compute->invoked_flag |= Compute::INVOKED_VECTOR;
    }
};

TEST_F(FixUVTTest, ComputeDependencyIsRefreshedAfterDrift)
{
    setup_quadratic_system();
    command("velocity all set 1.0 0.0 0.0");
    command("compute position all reduce sum x y");
    command("variable response equal c_position[1]");
    (*lmp->modify->fix_map)["early/compute"] =
        [](LAMMPS *lmp, int narg, char **arg) -> Fix * {
            return new FixEarlyCompute(lmp, narg, arg);
        };
    command("fix early all early/compute");
    command("fix cp all uvt temp 1 1 100000 mu 2 2 0.5 ne 1.8 dedn v_response");
    command("run 1 post no");
    double local = 0.0, total = 0.0;
    for (int i = 0; i < lmp->atom->nlocal; ++i) local += lmp->atom->x[i][0];
    MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_NEAR(fix_value("cp", 14), total, 1.0e-10);
}

TEST_F(FixUVTTest, RespaElectronicKicksUseOutermostTimestep)
{
    setup_quadratic_system();
    command("run_style respa 2 4");
    command("variable response equal 4.0");
    command("fix cp all uvt temp 1 1 1.0e10 mu 2 2 0.5 ne 1.8 dedn v_response");
    command("run 0 post no");
    int dim = 0;
    auto *cp = lmp->modify->get_fix_by_id("cp");
    const double mass = *static_cast<double *>(cp->extract("ne_mass", dim));
    command("run 1 post no");
    EXPECT_NEAR(fix_value("cp", 13), -2.0 * 0.005 / mass, 1.0e-12);
    EXPECT_NEAR(fix_value("cp", 12), 1.8 - 0.005 * 0.005 / mass, 1.0e-12);
    EXPECT_DOUBLE_EQ(fix_value("cp", 14), 4.0);
}

TEST_F(FixUVTTest, QuadraticToyPhysicsAveragesConverge)
{
    setup_quadratic_system();
    setup_quadratic_fix();
    command("fix avg all ave/time 1 10000 10000 f_cp[13] f_cp[14] f_cp[15]");
    command("run 10000 post no");

    EXPECT_NEAR(fix_value("avg", 0), 1.4, 1.0e-1);
    EXPECT_NEAR(fix_value("avg", 1), 0.0, 1.0e-1);
    EXPECT_NEAR(fix_value("avg", 2), 2.0, 1.0e-1);
}

TEST_F(FixUVTTest, RestartRestoresElectronState)
{
    double ne_before = 0.0;
    double nedot_before = 0.0;
    double dedn_before = 0.0;
    double energy_before = 0.0;
    setup_quadratic_system();
    setup_quadratic_fix();
    command("run 1000 post no");
    ne_before = fix_value("cp", 12);
    nedot_before = fix_value("cp", 13);
    dedn_before = fix_value("cp", 14);
    energy_before = fix_value("cp", 15);
    command("write_restart uvt.restart");
    command("clear");
    command("read_restart uvt.restart");
    command("variable k_quad equal 5.0");
    command("variable N0_quad equal 1.0");
    command("variable dEdN equal v_k_quad*(f_cp[13]-v_N0_quad)");
    command("fix cp all uvt temp 1.0 1.0 0.5 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
    command("run 0 post no");
    platform::unlink("uvt.restart");

    EXPECT_NEAR(fix_value("cp", 12), ne_before, 1.0e-10);
    EXPECT_NEAR(fix_value("cp", 13), nedot_before, 1.0e-10);
    EXPECT_NEAR(fix_value("cp", 14), dedn_before, 1.0e-10);
    EXPECT_NEAR(fix_value("cp", 15), energy_before, 1.0e-10);
}

TEST_F(FixUVTTest, ExtractCurrentMuReportsCurrentDEDN)
{
    setup_quadratic_system();
    command("variable k_quad equal 5.0");
    command("variable N0_quad equal 1.0");
    command("variable dEdN equal v_k_quad*(f_cp[13]-v_N0_quad)");
    command("fix cp all uvt temp 1.0 1.0 0.5 mu 2.0 3.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
    command("run 0 post no");

    auto *fix = lmp->modify->get_fix_by_id("cp");
    ASSERT_NE(fix, nullptr);

    int dim = 0;
    auto *u_current = static_cast<double *>(fix->extract("u_current", dim));
    ASSERT_NE(u_current, nullptr);
    EXPECT_EQ(dim, 0);
    auto *dedn = static_cast<double *>(fix->extract("dedn", dim));
    EXPECT_EQ(u_current, dedn);
    EXPECT_EQ(dim, 1);

    EXPECT_NEAR(*u_current, fix_value("cp", 14), 1.0e-12);
    EXPECT_NE(*u_current, fix_value("cp", 15));
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
