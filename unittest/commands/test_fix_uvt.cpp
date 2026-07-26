/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#include "lammps.h"
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
