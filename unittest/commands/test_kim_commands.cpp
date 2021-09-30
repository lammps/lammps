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
#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "lmppython.h"
#include "modify.h"
#include "output.h"
#include "utils.h"
#include "variable.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::MatchesRegex;
using ::testing::StrEq;

class KimCommandsTest : public LAMMPSTest {
protected:
    Variable *variable;

    void SetUp() override
    {
        testbinary = "KimCommandsTest";
        LAMMPSTest::SetUp();
        variable = lmp->input->variable;
    }
};

TEST_F(KimCommandsTest, kim)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal kim command.*", command("kim"););
    TEST_FAILURE(".*ERROR: Unknown kim subcommand.*", command("kim unknown"););
    TEST_FAILURE(".*kim_init.*has been renamed to.*", command("kim_init"););
    TEST_FAILURE(".*kim_interactions.*has been renamed to.*", command("kim_interactions"););
    TEST_FAILURE(".*kim_param.*has been renamed to.*", command("kim_param"););
    TEST_FAILURE(".*kim_property.*has been renamed to.*", command("kim_property"););
    TEST_FAILURE(".*kim_query.*has been renamed to.*", command("kim_query"););
}

TEST_F(KimCommandsTest, kim_init)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal 'kim init' command.*", command("kim init"););
    TEST_FAILURE(".*ERROR: Illegal 'kim init' command.*",
                 command("kim init LennardJones_Ar real si"););
    TEST_FAILURE(".*ERROR: LAMMPS unit_style lj not supported by KIM models.*",
                 command("kim init LennardJones_Ar lj"););
    TEST_FAILURE(".*ERROR: LAMMPS unit_style micro not supported by KIM models.*",
                 command("kim init LennardJones_Ar micro"););
    TEST_FAILURE(".*ERROR: LAMMPS unit_style nano not supported by KIM models.*",
                 command("kim init LennardJones_Ar nano"););
    TEST_FAILURE(".*ERROR: Unknown unit_style.*", command("kim init LennardJones_Ar new_style"););
    TEST_FAILURE(".*ERROR: KIM Model name not found.*", command("kim init Unknown_Model real"););
    TEST_FAILURE(".*ERROR: Incompatible units for KIM Simulator Model, required units = metal.*",
                 command("kim init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu real"););
    // TEST_FAILURE(".*ERROR: KIM Model does not support the requested unit system.*",
    //              command("kim init ex_model_Ar_P_Morse real"););

    BEGIN_HIDE_OUTPUT();
    command("kim init LennardJones_Ar real");
    END_HIDE_OUTPUT();

    int ifix = lmp->modify->find_fix("KIM_MODEL_STORE");
    ASSERT_GE(ifix, 0);
}

TEST_F(KimCommandsTest, kim_interactions)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal 'kim interactions' command.*", command("kim interactions"););

    BEGIN_HIDE_OUTPUT();
    command("kim init LennardJones_Ar real");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Must use 'kim interactions' command "
                 "after simulation box is defined.*",
                 command("kim interactions Ar"););

    BEGIN_HIDE_OUTPUT();
    command("kim init LennardJones_Ar real");
    command("lattice fcc 4.4300");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Illegal 'kim interactions' command.*",
                 command("kim interactions Ar Ar"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("lattice fcc 4.4300");
    command("region box block 0 20 0 20 0 20");
    command("create_box 4 box");
    command("create_atoms 4 box");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Illegal 'kim interactions' command.*",
                 command("kim interactions Ar Ar"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("lattice fcc 4.4300");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Must use 'kim init' before 'kim interactions'.*",
                 command("kim interactions Ar"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init LennardJones_Ar real");
    command("lattice fcc 4.4300");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: fixed_types cannot be used with a KIM Portable Model.*",
                 command("kim interactions fixed_types"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("pair_style kim LennardJones_Ar");
    command("region box block 0 1 0 1 0 1");
    command("create_box 4 box");
    command("pair_coeff * * Ar Ar Ar Ar");
    END_HIDE_OUTPUT();

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu metal");
    command("lattice fcc 4.920");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Species 'Ar' is not supported by this KIM Simulator Model.*",
                 command("kim interactions Ar"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu metal");
    command("lattice fcc 4.08");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("kim interactions Au");
    END_HIDE_OUTPUT();

    // ASSERT_EQ(lmp->output->var_kim_periodic, 1);
    // TEST_FAILURE(".*ERROR: Incompatible units for KIM Simulator Model.*",
    //              command("kim interactions Au"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init LennardJones_Ar real");
    command("lattice fcc 4.4300");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("kim interactions Ar");
    command("mass 1 39.95");
    END_HIDE_OUTPUT();

    int ifix = lmp->modify->find_fix("KIM_MODEL_STORE");
    ASSERT_GE(ifix, 0);

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init LennardJones_Ar real");
    command("lattice fcc 4.4300");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("kim interactions Ar");
    command("mass 1 39.95");
    command("run 1");
    command("kim interactions Ar");
    command("run 1");
    END_HIDE_OUTPUT();
}

TEST_F(KimCommandsTest, kim_param)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal 'kim param' command.*", command("kim param"););
    TEST_FAILURE(".*ERROR: Incorrect arguments in 'kim param' command.\n"
                 "'kim param get/set' is mandatory.*",
                 command("kim param unknown shift 1 shift"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu metal");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: 'kim param' can only be used with a KIM Portable Model.*",
                 command("kim param get shift 1 shift"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init LennardJones612_UniversalShifted__MO_959249795837_003 real");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Illegal 'kim param get' command.\nTo get the new "
                 "parameter values, pair style must be assigned.\nMust use 'kim"
                 " interactions' or 'pair_style kim' before 'kim param get'.*",
                 command("kim param get shift 1 shift"););

    TEST_FAILURE(".*ERROR: Illegal 'kim param set' command.\nTo set the new "
                 "parameter values, pair style must be assigned.\nMust use 'kim"
                 " interactions' or 'pair_style kim' before 'kim param set'.*",
                 command("kim param set shift 1 2"););

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init LennardJones612_UniversalShifted__MO_959249795837_003 real");
    command("lattice fcc 4.4300");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("kim interactions Ar");
    command("mass 1 39.95");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Illegal index '0' for "
                 "'shift' parameter with the extent of '1'.*",
                 command("kim param get shift 0 shift"););
    TEST_FAILURE(".*ERROR: Illegal index '2' for "
                 "'shift' parameter with the extent of '1'.*",
                 command("kim param get shift 2 shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range.\nExpected integer "
                 "parameter\\(s\\) instead of '1.' in index_range.*",
                 command("kim param get shift 1. shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range '1-2' for 'shift' "
                 "parameter with the extent of '1'.*",
                 command("kim param get shift 1:2 shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range.\nExpected integer "
                 "parameter\\(s\\) instead of '1-2' in index_range.*",
                 command("kim param get shift 1-2 shift"););
    TEST_FAILURE(".*ERROR: Wrong number of arguments in 'kim param "
                 "get' command.\nThe LAMMPS '3' variable names or "
                 "'s1 split' is mandatory.*",
                 command("kim param get sigmas 1:3 s1 s2"););
    TEST_FAILURE(".*ERROR: Wrong argument in 'kim param get' command.\nThis "
                 "Model does not have the requested 'unknown' parameter.*",
                 command("kim param get unknown 1 unknown"););

    BEGIN_HIDE_OUTPUT();
    command("kim param get shift 1 shift");
    END_HIDE_OUTPUT();

    ASSERT_FALSE(variable->find("shift") == -1);
    ASSERT_THAT(variable->retrieve("shift"), StrEq("1"));

    TEST_FAILURE(".*ERROR: Illegal index '2' for "
                 "'shift' parameter with the extent of '1'.*",
                 command("kim param set shift 2 2"););
    TEST_FAILURE(".*ERROR: Illegal index_range.\nExpected integer "
                 "parameter\\(s\\) instead of '1.' in index_range.*",
                 command("kim param set shift 1. shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range '1-2' for "
                 "'shift' parameter with the extent of '1'.*",
                 command("kim param set shift 1:2 2"););
    TEST_FAILURE(".*ERROR: Wrong number of variable values for pair coefficients.*",
                 command("kim param set sigmas 1:3 0.5523570 0.4989030"););
    TEST_FAILURE(".*ERROR: Wrong argument for pair coefficients.\nThis "
                 "Model does not have the requested '0.4989030' parameter.*",
                 command("kim param set sigmas 1:1 0.5523570 0.4989030"););

    BEGIN_HIDE_OUTPUT();
    command("variable new_shift equal 2");
    command("kim param set shift 1 ${new_shift}");
    command("kim param get shift 1 shift");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("shift"), StrEq("2"));

    TEST_FAILURE(".*ERROR: Illegal variable name in 'kim param get'.*",
                 command("kim param get cutoffs 1:3 list"););
    TEST_FAILURE(".*ERROR: Illegal variable name in 'kim param get'.*",
                 command("kim param get cutoffs 1:3 cutoffs_1 cutoffs_2 list"););
    TEST_FAILURE(".*ERROR: Illegal variable name in 'kim param get'.*",
                 command("kim param get cutoffs 1:3 split"););
    TEST_FAILURE(".*ERROR: Illegal variable name in 'kim param get'.*",
                 command("kim param get cutoffs 1:3 cutoffs_1 cutoffs_2 split"););
    TEST_FAILURE(".*ERROR: Illegal variable name in 'kim param get'.*",
                 command("kim param get cutoffs 1:3 explicit"););
    TEST_FAILURE(".*ERROR: Illegal variable name in 'kim param get'.*",
                 command("kim param get cutoffs 1:3 cutoffs_1 cutoffs_2 explicit"););
    TEST_FAILURE(".*ERROR: Wrong number of arguments in 'kim param get' "
                 "command.\nThe LAMMPS '3' variable names or 'cutoffs "
                 "split/list' is mandatory.*",
                 command("kim param get cutoffs 1:3 cutoffs"););
    TEST_FAILURE(".*ERROR: Wrong number of arguments in 'kim param get' "
                 "command.\nThe LAMMPS '3' variable names or 'cutoffs_1 "
                 "split' is mandatory.*",
                 command("kim param get cutoffs 1:3 cutoffs_1 cutoffs_2"););

    BEGIN_HIDE_OUTPUT();
    command("kim param get cutoffs 1:3 cutoffs_1 cutoffs_2 cutoffs_3");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("cutoffs_1"), StrEq("2.20943"));
    ASSERT_THAT(variable->retrieve("cutoffs_2"), StrEq("2.10252"));
    ASSERT_THAT(variable->retrieve("cutoffs_3"), StrEq("5.666115"));

    BEGIN_HIDE_OUTPUT();
    command("kim param get cutoffs 1:3 cutoffs_1 cutoffs_2 cutoffs_3 explicit");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("cutoffs_1"), StrEq("2.20943"));
    ASSERT_THAT(variable->retrieve("cutoffs_2"), StrEq("2.10252"));
    ASSERT_THAT(variable->retrieve("cutoffs_3"), StrEq("5.666115"));

    BEGIN_HIDE_OUTPUT();
    command("kim param get cutoffs 1:3 cutoffs split");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("cutoffs_1"), StrEq("2.20943"));
    ASSERT_THAT(variable->retrieve("cutoffs_2"), StrEq("2.10252"));
    ASSERT_THAT(variable->retrieve("cutoffs_3"), StrEq("5.666115"));

    BEGIN_HIDE_OUTPUT();
    command("kim param get cutoffs 1:3 cutoffs list");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("cutoffs"), StrEq("2.20943 2.10252 5.666115"));

    BEGIN_HIDE_OUTPUT();
    command("kim param set cutoffs 1 2.21 cutoffs 2 2.11");
    command("kim param get cutoffs 1:2 cutoffs list");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("cutoffs"), StrEq("2.21 2.11"));

    BEGIN_HIDE_OUTPUT();
    command("kim param set cutoffs 1:3 2.3 2.2 5.7");
    command("kim param get cutoffs 1:3 cutoffs list");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("cutoffs"), StrEq("2.3 2.2 5.7"));

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    command("lattice fcc 4.4300");
    command("region box block 0 10 0 10 0 10");
    command("create_box 1 box");
    command("create_atoms 1 box");
    command("mass 1 39.95");
    command("pair_style kim LennardJones612_UniversalShifted__MO_959249795837_003");
    command("pair_coeff * * Ar");
    END_HIDE_OUTPUT();

    BEGIN_HIDE_OUTPUT();
    command("kim param get shift 1 shift");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("shift"), StrEq("1"));

    BEGIN_HIDE_OUTPUT();
    command("variable new_shift equal 2");
    command("kim param set shift 1 ${new_shift}");
    command("kim param get shift 1 shift");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("shift"), StrEq("2"));
}

TEST_F(KimCommandsTest, kim_property)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();
    if (!LAMMPS::is_installed_pkg("PYTHON")) GTEST_SKIP();

    if (!lmp->python->has_minimum_version(3, 6)) {
        TEST_FAILURE(".*ERROR: Invalid Python version.\n"
                     "The kim-property Python package requires Python "
                     "3 >= 3.6 support.*",
                     command("kim property"););
    } else {
        TEST_FAILURE(".*ERROR: Invalid 'kim property' command.*", command("kim property"););
        TEST_FAILURE(".*ERROR: Invalid 'kim property' command.*", command("kim property create"););
        TEST_FAILURE(".*ERROR: Incorrect arguments in 'kim property' command."
                     "\n'kim property create/destroy/modify/remove/dump' "
                     "is mandatory.*",
                     command("kim property unknown 1 atomic-mass"););
    }
#if defined(KIM_EXTRA_UNITTESTS)
    TEST_FAILURE(".*ERROR: Invalid 'kim property create' command.*",
                 command("kim property create 1"););
    TEST_FAILURE(".*ERROR: Invalid 'kim property destroy' command.*",
                 command("kim property destroy 1 cohesive-potential-energy-cubic-crystal"););
    TEST_FAILURE(".*ERROR: Invalid 'kim property modify' command.*",
                 command("kim property modify 1 key short-name"););
    TEST_FAILURE(".*ERROR: There is no property instance to modify the content.*",
                 command("kim property modify 1 key short-name source-value 1 fcc"););
    TEST_FAILURE(".*ERROR: Invalid 'kim property remove' command.*",
                 command("kim property remove 1 key"););
    TEST_FAILURE(".*ERROR: There is no property instance to remove the content.*",
                 command("kim property remove 1 key short-name"););
    TEST_FAILURE(".*ERROR: There is no property instance to dump the content.*",
                 command("kim property dump results.edn"););
    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init LennardJones612_UniversalShifted__MO_959249795837_003 real");
    command("kim property create 1 cohesive-potential-energy-cubic-crystal");
    command("kim property modify 1 key short-name source-value 1 fcc");
    command("kim property destroy 1");
    END_HIDE_OUTPUT();
#endif
}

TEST_F(KimCommandsTest, kim_query)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.*", command("kim query"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe keyword 'split' "
                 "must be followed by the name of the query function.*",
                 command("kim query a0 split"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe keyword 'list' "
                 "must be followed by the name of the query function.*",
                 command("kim query a0 list"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe keyword 'index' "
                 "must be followed by the name of the query function.*",
                 command("kim query a0 index"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe 'list' keyword "
                 "can not be used after 'split'.*",
                 command("kim query a0 split list"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe 'index' keyword "
                 "can not be used after 'split'.*",
                 command("kim query a0 split index"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe 'split' keyword "
                 "can not be used after 'list'.*",
                 command("kim query a0 list split"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe 'index' keyword "
                 "can not be used after 'list'.*",
                 command("kim query a0 list index"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe 'list' keyword "
                 "can not be used after 'index'.*",
                 command("kim query a0 index list"););
    TEST_FAILURE(".*ERROR: Illegal 'kim query' command.\nThe 'split' keyword "
                 "can not be used after 'index'.*",
                 command("kim query a0 index split"););
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `crystal` "
                 "to 'kim query' is wrong. The query format is the "
                 "keyword=\\[value\\], where value is always an array of one "
                 "or more comma-separated items.*",
                 command("kim query a0 get_lattice_constant_cubic "
                         "crystal"););
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `"
                 "crystal=fcc` to 'kim query' is wrong. The query format is the "
                 "keyword=\\[value\\], where value is always an array of one "
                 "or more comma-separated items.*",
                 command("kim query a0 get_lattice_constant_cubic "
                         "crystal=fcc"););
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `"
                 "crystal=\\[fcc` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command("kim query a0 get_lattice_constant_cubic "
                         "crystal=[fcc"););
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `"
                 "crystal=fcc\\]` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command("kim query a0 get_lattice_constant_cubic "
                         "crystal=fcc]"););

    std::string squery = "kim query a0 get_lattice_constant_cubic "
                         "crystal=[\"fcc\"] species=\"Al\",\"Ni\" units=[\"angstrom\"]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "\"Al\",\"Ni\"` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command(squery););

    squery = "kim query a0 get_lattice_constant_cubic crystal=[fcc] species=Al,Ni units=[angstrom]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "Al,Ni` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command(squery););

    squery =
        "kim query a0 get_lattice_constant_cubic crystal=[fcc] species=Al,Ni, units=[angstrom]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "Al,Ni,` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command(squery););

    squery =
        "kim query a0 get_lattice_constant_cubic crystal=[fcc] species=[Al,Ni, units=[angstrom]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "\\[Al,Ni,` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command(squery););

    squery =
        "kim query a0 get_lattice_constant_cubic crystal=[fcc] species=Al,Ni], units=[angstrom]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "Al,Ni\\],` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command(squery););

    squery = "kim query a0 get_lattice_constant_cubic crystal=[fcc] species=Al,\"Ni\"], "
             "units=[angstrom]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "Al,\"Ni\"\\],` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command(squery););

    squery = "kim query a0 get_lattice_constant_cubic crystal=[fcc] species=\"Al\",Ni], "
             "units=[angstrom]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "\"Al\",Ni\\],` to 'kim query' is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*",
                 command(squery););

    squery = "kim query a0 get_lattice_constant_cubic crystal=[\"fcc\"] species=[\"Al\"]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nMust use 'kim init' before "
                 "'kim query' or must provide the model name after query "
                 "function with the format of 'model=\\[model_name\\]'.*",
                 command(squery););

    squery = "kim query a0 get_lattice_constant_cubic crystal=[fcc] species=[Al]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nMust use 'kim init' before "
                 "'kim query' or must provide the model name after query "
                 "function with the format of 'model=\\[model_name\\]'.*",
                 command(squery););

    squery = "kim query a0 get_lattice_constant_cubic crystal=[\"fcc\"] species=[Al]";
    TEST_FAILURE(".*ERROR: Illegal query format.\nMust use 'kim init' before "
                 "'kim query' or must provide the model name after query "
                 "function with the format of 'model=\\[model_name\\]'.*",
                 command(squery););

#if defined(KIM_EXTRA_UNITTESTS)
    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim query latconst_1 get_lattice_constant_cubic "
            "crystal=[fcc] species=[Al] units=[angstrom] "
            "model=[EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005]");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("latconst_1"), StrEq("4.032082033157349"));

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal");
    command("kim query latconst_1 get_lattice_constant_cubic crystal=[fcc] species=[Al] "
            "units=[angstrom]");

    command("kim query latconst_2 get_lattice_constant_cubic crystal=[fcc] species=[Al] "
            "units=[angstrom] "
            "model=[LennardJones612_UniversalShifted__MO_959249795837_003]");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("latconst_1"), StrEq("4.032082033157349"));
    ASSERT_THAT(variable->retrieve("latconst_2"), StrEq("3.328125931322575"));

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init EAM_Dynamo_MendelevAckland_2007v3_Zr__MO_004835508849_000 metal");

    command("kim query latconst split get_lattice_constant_hexagonal crystal=[hcp] species=[Zr] "
            "units=[angstrom]");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("latconst_1"), StrEq("3.234055244384789"));
    ASSERT_THAT(variable->retrieve("latconst_2"), StrEq("5.167650199630013"));

    BEGIN_HIDE_OUTPUT();
    command("clear");

    command("kim query latconst index get_lattice_constant_hexagonal "
            "crystal=[hcp] species=[Zr] units=[angstrom] "
            "model=[EAM_Dynamo_MendelevAckland_2007v3_Zr__MO_004835508849_000]");
    END_HIDE_OUTPUT();
    ASSERT_THAT(variable->retrieve("latconst"), StrEq("3.234055244384789"));

    BEGIN_HIDE_OUTPUT();
    command("variable latconst delete");
    command("clear");
    command("kim init EAM_Dynamo_MendelevAckland_2007v3_Zr__MO_004835508849_000 metal");

    command("kim query latconst list get_lattice_constant_hexagonal crystal=[hcp] species=[Zr] "
            "units=[angstrom]");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("latconst"), StrEq("3.234055244384789 5.167650199630013"));

    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("kim init EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal");

    command("kim query alpha get_linear_thermal_expansion_coefficient_cubic "
            "crystal=[fcc] species=[Al] units=[1/K] temperature=[293.15] "
            "temperature_units=[K]");
    END_HIDE_OUTPUT();

    ASSERT_THAT(variable->retrieve("alpha"), StrEq("1.654960564704273e-05"));

    BEGIN_HIDE_OUTPUT();
    command("clear");

    command("kim query model_list list get_available_models species=[Al]");
    END_HIDE_OUTPUT();

    std::string model_list = variable->retrieve("model_list");
    auto n = model_list.find("EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005");
    ASSERT_TRUE(n != std::string::npos);

    BEGIN_HIDE_OUTPUT();
    command("clear");

    command("kim query model_name index get_available_models species=[Al]");
    command("variable model_name delete");
    END_HIDE_OUTPUT();

    BEGIN_HIDE_OUTPUT();
    command("clear");

    command("kim query model_name index get_available_models "
            "species=[Al] potential_type=[eam,meam]");
    command("variable model_name delete");

    command("kim query model_name index get_available_models "
            "species=[Al] potential_type=[\"eam\",\"meam\"]");
    command("variable model_name delete");

    command("kim query model_name index get_available_models "
            "species=[Al] potential_type=[eam,\"meam\"]");
    command("variable model_name delete");

    command("kim query model_name index get_available_models "
            "species=[Al] potential_type=[\"eam\",meam]");
    command("variable model_name delete");
    END_HIDE_OUTPUT();
#endif
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (Info::get_mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

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
    return rv;
}
