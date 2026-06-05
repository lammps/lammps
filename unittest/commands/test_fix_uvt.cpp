/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#include "lammps.h"

#include "platform.h"

#include "../testing/core.h"
#include "gtest/gtest.h"

#include <mpi.h>

bool verbose = false;

namespace LAMMPS_NS {

class FixUVTTest : public LAMMPSTest {
protected:
    void setup_quadratic_system()
    {
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
        command("variable mu_target equal 2.0");
        command("variable k_quad equal 5.0");
        command("variable N0_quad equal 1.0");
        command("variable dEdN equal v_k_quad*(v_Ne-v_N0_quad)");
        command("fix cp all uvt temp 1.0 1.0 0.5 mu 2.0 2.0 0.5 ne 1.8 ne_velocity 0.0 dedn v_dEdN");
        command("variable Ne equal f_cp[13]");
        command("variable NeDot equal f_cp[14]");
        command("variable dEdNfix equal f_cp[15]");
        command("fix avg all ave/time 1 1000 1000 v_Ne v_NeDot v_dEdNfix mode vector");
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
    BEGIN_HIDE_OUTPUT();
    setup_quadratic_system();
    setup_quadratic_fix();
    command("run 1000 post no");
    END_HIDE_OUTPUT();

    EXPECT_NEAR(fix_value("avg", 0), 1.4, 5.0e-2);
    EXPECT_NEAR(fix_value("avg", 1), 0.0, 5.0e-2);
    EXPECT_NEAR(fix_value("avg", 2), 2.0, 5.0e-2);
}

TEST_F(FixUVTTest, RestartRestoresElectronState)
{
    BEGIN_HIDE_OUTPUT();
    setup_quadratic_system();
    setup_quadratic_fix();
    command("run 1000 post no");
    double ne_before = fix_value("cp", 12);
    double nedot_before = fix_value("cp", 13);
    double dedn_before = fix_value("cp", 14);
    double energy_before = fix_value("cp", 15);
    command("write_restart uvt.restart");
    command("clear");
    command("read_restart uvt.restart");
    platform::unlink("uvt.restart");
    END_HIDE_OUTPUT();

    EXPECT_NEAR(fix_value("cp", 12), ne_before, 1.0e-12);
    EXPECT_NEAR(fix_value("cp", 13), nedot_before, 1.0e-12);
    EXPECT_NEAR(fix_value("cp", 14), dedn_before, 1.0e-12);
    EXPECT_NEAR(fix_value("cp", 15), energy_before, 1.0e-12);
}

}    // namespace LAMMPS_NS
