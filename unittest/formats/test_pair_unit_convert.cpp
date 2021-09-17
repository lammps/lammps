/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/core.h"
#include "atom.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "output.h"
#include "pair.h"
#include "thermo.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <vector>

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

class PairUnitConvertTest : public LAMMPSTest {
protected:
    double fold[4][3];

    void SetUp() override
    {
        testbinary = "PairUnitConvertTest";
        LAMMPSTest::SetUp();
        ASSERT_NE(lmp, nullptr);

        BEGIN_HIDE_OUTPUT();
        command("units metal");
        command("dimension 3");
        command("region box block -4 4 -4 4 -4 4");
        command("create_box 2 box");
        command("create_atoms 1 single -1.1  1.2  0.0 units box");
        command("create_atoms 1 single -1.2 -1.1  0.0 units box");
        command("create_atoms 2 single  0.9  1.0  0.0 units box");
        command("create_atoms 2 single  1.0 -0.9  0.0 units box");
        command("pair_style zero 4.0");
        command("pair_coeff * *");
        command("mass * 1.0");
        command("write_data test_pair_unit_convert.data nocoeff");
        command("clear");
        END_HIDE_OUTPUT();
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        remove("test_pair_unit_convert.data");
    }
};

TEST_F(PairUnitConvertTest, zero)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "zero")) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style zero 6.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style zero 6.0");
    command("pair_coeff * *");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style lj/cut 6.0");
    command("pair_coeff * * 0.01014286346782117 2.0");
    remove("test.table.metal");
    command("pair_write 1 1 1000 r 0.1 6.0 test.table.metal lj_1_1");
    command("pair_write 1 2 1000 r 0.1 6.0 test.table.metal lj_1_2");
    command("pair_write 2 2 1000 r 0.1 6.0 test.table.metal lj_2_2");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style lj/cut 6.0");
    command("pair_coeff * * 0.2339 2.0");
    remove("test.table.real");
    command("pair_write 1 1 1000 r 0.1 6.0 test.table.real lj_1_1");
    command("pair_write 1 2 1000 r 0.1 6.0 test.table.real lj_1_2");
    command("pair_write 2 2 1000 r 0.1 6.0 test.table.real lj_2_2");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam");
    command("pair_coeff * * Cu_u3.eam");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam");
    command("pair_coeff * * Cu_u3.eam");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam/alloy");
    command("pair_coeff * * AlCu.eam.alloy Al Cu");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam/alloy");
    command("pair_coeff * * AlCu.eam.alloy Al Cu");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam/fs");
    command("pair_coeff * * FeP_mm.eam.fs Fe P");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam/fs");
    command("pair_coeff * * FeP_mm.eam.fs Fe P");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam/cd");
    command("pair_coeff * * FeCr.cdeam Cr Fe");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eam/cd");
    command("pair_coeff * * FeCr.cdeam Cr Fe");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eim");
    command("pair_coeff * * Na Cl ffield.eim Na Cl");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style eim");
    command("pair_coeff * * Na Cl ffield.eim Na Cl");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style gw");
    command("pair_coeff * * SiC.gw Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style gw");
    command("pair_coeff * * SiC.gw Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style gw/zbl");
    command("pair_coeff * * SiC.gw.zbl Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style gw/zbl");
    command("pair_coeff * * SiC.gw.zbl Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style nb3b/harmonic");
    command("pair_coeff * * MOH.nb3b.harmonic M O");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style nb3b/harmonic");
    command("pair_coeff * * MOH.nb3b.harmonic M O");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style sw");
    command("pair_coeff * * GaN.sw Ga N");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style sw");
    command("pair_coeff * * GaN.sw Ga N");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style table linear 1000");
    command("pair_coeff 1 1 test.table.metal lj_1_1");
    command("pair_coeff 1 2 test.table.metal lj_1_2");
    command("pair_coeff 2 2 test.table.metal lj_2_2");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style table linear 1000");
    command("pair_coeff 1 1 test.table.metal lj_1_1");
    command("pair_coeff 1 2 test.table.metal lj_1_2");
    command("pair_coeff 2 2 test.table.metal lj_2_2");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style table linear 1000");
    command("pair_coeff 1 1 test.table.real lj_1_1");
    command("pair_coeff 1 2 test.table.real lj_1_2");
    command("pair_coeff 2 2 test.table.real lj_2_2");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style table linear 1000");
    command("pair_coeff 1 1 test.table.real lj_1_1");
    command("pair_coeff 1 2 test.table.real lj_1_2");
    command("pair_coeff 2 2 test.table.real lj_2_2");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    double pnew;
    lmp->output->thermo->evaluate_keyword("press", &pnew);
    EXPECT_NEAR(pold, 1.0 / p_convert * pnew, fabs(pnew * rel_error));
    double enew = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    EXPECT_NEAR(1.0 / ev_convert * eold, enew, fabs(enew * rel_error));

    f = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(1.0 / ev_convert * fold[i][j], f[i][j], fabs(f[i][j] * rel_error));
}

TEST_F(PairUnitConvertTest, tersoff)
{
    // check if the prerequisite pair style is available
    if (!info->has_style("pair", "tersoff")) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff");
    command("pair_coeff * * SiC.tersoff Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff");
    command("pair_coeff * * SiC.tersoff Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/mod");
    command("pair_coeff * * Si.tersoff.mod Si Si");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/mod");
    command("pair_coeff * * Si.tersoff.mod Si Si");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/mod/c");
    command("pair_coeff * * Si.tersoff.modc Si Si");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/mod/c");
    command("pair_coeff * * Si.tersoff.modc Si Si");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/table");
    command("pair_coeff * * SiC.tersoff Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/table");
    command("pair_coeff * * SiC.tersoff Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/zbl");
    command("pair_coeff * * SiC.tersoff.zbl Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/zbl");
    command("pair_coeff * * SiC.tersoff.zbl Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("package omp 4");
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/zbl/omp");
    command("pair_coeff * * SiC.tersoff.zbl Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("package omp 4");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style tersoff/zbl/omp");
    command("pair_coeff * * SiC.tersoff.zbl Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("read_data test_pair_unit_convert.data");
    command("pair_style vashishta");
    command("pair_coeff * * SiC.vashishta Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    // copy pressure, energy, and force from first step
    double pold;
    lmp->output->thermo->evaluate_keyword("press", &pold);
    double eold = lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
    double **f  = lmp->atom->f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            fold[i][j] = f[i][j];

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("read_data test_pair_unit_convert.data");
    command("pair_style vashishta");
    command("pair_coeff * * SiC.vashishta Si C");
    command("run 0 post no");
    END_HIDE_OUTPUT();

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
