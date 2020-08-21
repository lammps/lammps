/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "output.h"
#include "pair.h"
#include "thermo.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::Eq;

// eV to kcal/mol conversion constant (CODATA 2018)
const double ev_convert = utils::get_conversion_factor(utils::ENERGY, utils::METAL2REAL);
// 1atm in bar
const double p_convert = 1.01325;
// relative error for comparing numbers
// cannot use smaller value due to lack of consistency
// of data in update.cpp. could be 1.0e-12
const double rel_error = 5.0e-7;

class PairUnitConvertTest : public ::testing::Test {
protected:
    LAMMPS *lmp;
    Info *info;
    double fold[4][3];

    void SetUp() override
    {
        const char *args[] = {"PairUnitConvertTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        ASSERT_NE(lmp, nullptr);
        if (!verbose) ::testing::internal::CaptureStdout();
        info = new Info(lmp);
        lmp->input->one("units metal");
        lmp->input->one("dimension 3");
        lmp->input->one("region box block -4 4 -4 4 -4 4");
        lmp->input->one("create_box 2 box");
        lmp->input->one("create_atoms 1 single -1.1  1.2  0.0 units box");
        lmp->input->one("create_atoms 1 single -1.2 -1.1  0.0 units box");
        lmp->input->one("create_atoms 2 single  0.9  1.0  0.0 units box");
        lmp->input->one("create_atoms 2 single  1.0 -0.9  0.0 units box");
        lmp->input->one("pair_style zero 4.0");
        lmp->input->one("pair_coeff * *");
        lmp->input->one("mass * 1.0");
        lmp->input->one("write_data test_pair_unit_convert.data nocoeff");
        lmp->input->one("clear");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete info;
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
        remove("test_pair_unit_convert.data");
    }
};

TEST_F(PairUnitConvertTest, zero)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "zero")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style zero 6.0");
    lmp->input->one("pair_coeff * *");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style zero 6.0");
    lmp->input->one("pair_coeff * *");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, lj_cut)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "lj/cut")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style lj/cut 6.0");
    lmp->input->one("pair_coeff * * 0.01014286346782117 2.0");
    remove("test.table.metal");
    lmp->input->one("pair_write 1 1 1000 r 0.1 6.0 test.table.metal lj_1_1");
    lmp->input->one("pair_write 1 2 1000 r 0.1 6.0 test.table.metal lj_1_2");
    lmp->input->one("pair_write 2 2 1000 r 0.1 6.0 test.table.metal lj_2_2");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style lj/cut 6.0");
    lmp->input->one("pair_coeff * * 0.2339 2.0");
    remove("test.table.real");
    lmp->input->one("pair_write 1 1 1000 r 0.1 6.0 test.table.real lj_1_1");
    lmp->input->one("pair_write 1 2 1000 r 0.1 6.0 test.table.real lj_1_2");
    lmp->input->one("pair_write 2 2 1000 r 0.1 6.0 test.table.real lj_2_2");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, eam)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "eam")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam");
    lmp->input->one("pair_coeff * * Cu_u3.eam");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam");
    lmp->input->one("pair_coeff * * Cu_u3.eam");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, eam_alloy)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "eam/alloy")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam/alloy");
    lmp->input->one("pair_coeff * * AlCu.eam.alloy Al Cu");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam/alloy");
    lmp->input->one("pair_coeff * * AlCu.eam.alloy Al Cu");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, eam_fs)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "eam/fs")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam/fs");
    lmp->input->one("pair_coeff * * FeP_mm.eam.fs Fe P");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam/fs");
    lmp->input->one("pair_coeff * * FeP_mm.eam.fs Fe P");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, eam_cd)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "eam/cd")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam/cd");
    lmp->input->one("pair_coeff * * FeCr.cdeam Cr Fe");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eam/cd");
    lmp->input->one("pair_coeff * * FeCr.cdeam Cr Fe");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, eim)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "eim")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eim");
    lmp->input->one("pair_coeff * * Na Cl ffield.eim Na Cl");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style eim");
    lmp->input->one("pair_coeff * * Na Cl ffield.eim Na Cl");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, gw)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "gw")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style gw");
    lmp->input->one("pair_coeff * * SiC.gw Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style gw");
    lmp->input->one("pair_coeff * * SiC.gw Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, gw_zbl)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "gw/zbl")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style gw/zbl");
    lmp->input->one("pair_coeff * * SiC.gw.zbl Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style gw/zbl");
    lmp->input->one("pair_coeff * * SiC.gw.zbl Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, nb3b_harmonic)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "nb3b/harmonic")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style nb3b/harmonic");
    lmp->input->one("pair_coeff * * MOH.nb3b.harmonic M O");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style nb3b/harmonic");
    lmp->input->one("pair_coeff * * MOH.nb3b.harmonic M O");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, sw)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "sw")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style sw");
    lmp->input->one("pair_coeff * * GaN.sw Ga N");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style sw");
    lmp->input->one("pair_coeff * * GaN.sw Ga N");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, table_metal2real)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "table")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style table linear 1000");
    lmp->input->one("pair_coeff 1 1 test.table.metal lj_1_1");
    lmp->input->one("pair_coeff 1 2 test.table.metal lj_1_2");
    lmp->input->one("pair_coeff 2 2 test.table.metal lj_2_2");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style table linear 1000");
    lmp->input->one("pair_coeff 1 1 test.table.metal lj_1_1");
    lmp->input->one("pair_coeff 1 2 test.table.metal lj_1_2");
    lmp->input->one("pair_coeff 2 2 test.table.metal lj_2_2");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, table_real2metal)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "table")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style table linear 1000");
    lmp->input->one("pair_coeff 1 1 test.table.real lj_1_1");
    lmp->input->one("pair_coeff 1 2 test.table.real lj_1_2");
    lmp->input->one("pair_coeff 2 2 test.table.real lj_2_2");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style table linear 1000");
    lmp->input->one("pair_coeff 1 1 test.table.real lj_1_1");
    lmp->input->one("pair_coeff 1 2 test.table.real lj_1_2");
    lmp->input->one("pair_coeff 2 2 test.table.real lj_2_2");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, 1.0/p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(1.0/ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(1.0/ev_convert * fold[i][j], f[i][j],
                        fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, tersoff)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "tersoff")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff");
    lmp->input->one("pair_coeff * * SiC.tersoff Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff");
    lmp->input->one("pair_coeff * * SiC.tersoff Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, tersoff_mod)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "tersoff/mod")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/mod");
    lmp->input->one("pair_coeff * * Si.tersoff.mod Si Si");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/mod");
    lmp->input->one("pair_coeff * * Si.tersoff.mod Si Si");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, tersoff_mod_c)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "tersoff/mod/c")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/mod/c");
    lmp->input->one("pair_coeff * * Si.tersoff.modc Si Si");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/mod/c");
    lmp->input->one("pair_coeff * * Si.tersoff.modc Si Si");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, tersoff_table)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "tersoff/table")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/table");
    lmp->input->one("pair_coeff * * SiC.tersoff Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/table");
    lmp->input->one("pair_coeff * * SiC.tersoff Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, tersoff_zbl)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "tersoff/zbl")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/zbl");
    lmp->input->one("pair_coeff * * SiC.tersoff.zbl Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/zbl");
    lmp->input->one("pair_coeff * * SiC.tersoff.zbl Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, tersoff_zbl_omp)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "tersoff/zbl/omp")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("package omp 4");
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/zbl/omp");
    lmp->input->one("pair_coeff * * SiC.tersoff.zbl Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("package omp 4");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style tersoff/zbl/omp");
    lmp->input->one("pair_coeff * * SiC.tersoff.zbl Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, vashishta)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "vashishta")) GTEST_SKIP();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units metal");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style vashishta");
    lmp->input->one("pair_coeff * * SiC.vashishta Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("read_data test_pair_unit_convert.data");
    lmp->input->one("pair_style vashishta");
    lmp->input->one("pair_coeff * * SiC.vashishta Si C");
    lmp->input->one("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    remove("test.table.metal");
    remove("test.table.real");
    return rv;
}
