/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

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

#if defined(OMPI_MAJOR_VERSION)
const bool have_openmpi = true;
#else
const bool have_openmpi = false;
#endif

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::MatchesRegex;
using ::testing::StrEq;

#define TEST_FAILURE(errmsg, ...)                                 \
    if (Info::has_exceptions()) {                                 \
        ::testing::internal::CaptureStdout();                     \
        ASSERT_ANY_THROW({__VA_ARGS__});                          \
        auto mesg = ::testing::internal::GetCapturedStdout();     \
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (!have_openmpi) {                                      \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

class KimCommandsTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"KimCommandsTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(KimCommandsTest, kim_init)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal kim_init command.*", 
                 lmp->input->one("kim_init"););
    TEST_FAILURE(".*ERROR: Illegal kim_init command.*",
                 lmp->input->one("kim_init LennardJones_Ar real si"););
    TEST_FAILURE(".*ERROR: LAMMPS unit_style lj not supported by KIM models.*",
                 lmp->input->one("kim_init LennardJones_Ar lj"););
    TEST_FAILURE(".*ERROR: LAMMPS unit_style micro not supported by KIM models.*",
                 lmp->input->one("kim_init LennardJones_Ar micro"););
    TEST_FAILURE(".*ERROR: LAMMPS unit_style nano not supported by KIM models.*",
                 lmp->input->one("kim_init LennardJones_Ar nano"););
    TEST_FAILURE(".*ERROR: Unknown unit_style.*",
                 lmp->input->one("kim_init LennardJones_Ar new_style"););
    TEST_FAILURE(".*ERROR: KIM Model name not found.*",
                 lmp->input->one("kim_init Unknown_Model real"););
    TEST_FAILURE(".*ERROR: Incompatible units for KIM Simulator Model, required units = metal.*",
                 lmp->input->one("kim_init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu real"););
    // TEST_FAILURE(".*ERROR: KIM Model does not support the requested unit system.*",
    //              lmp->input->one("kim_init ex_model_Ar_P_Morse real"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("kim_init LennardJones_Ar real");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    int ifix = lmp->modify->find_fix("KIM_MODEL_STORE");
    ASSERT_GE(ifix, 0);
}

TEST_F(KimCommandsTest, kim_interactions)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal kim_interactions command.*",
                 lmp->input->one("kim_interactions"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("kim_init LennardJones_Ar real");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Must use 'kim_interactions' command "
                 "after simulation box is defined.*",
                 lmp->input->one("kim_interactions Ar"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("kim_init LennardJones_Ar real");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Illegal kim_interactions command.*",
                 lmp->input->one("kim_interactions Ar Ar"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 20 0 20 0 20");
    lmp->input->one("create_box 4 box");
    lmp->input->one("create_atoms 4 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Illegal kim_interactions command.*",
                 lmp->input->one("kim_interactions Ar Ar"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Must use 'kim_init' before 'kim_interactions'.*",
                 lmp->input->one("kim_interactions Ar"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init LennardJones_Ar real");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: fixed_types cannot be used with a KIM Portable Model.*",
                 lmp->input->one("kim_interactions fixed_types"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("units real");
    lmp->input->one("pair_style kim LennardJones_Ar");
    lmp->input->one("region box block 0 1 0 1 0 1");
    lmp->input->one("create_box 4 box");
    lmp->input->one("pair_coeff * * Ar Ar Ar Ar");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu metal");
    lmp->input->one("lattice fcc 4.920");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Species 'Ar' is not supported by this KIM Simulator Model.*",
                 lmp->input->one("kim_interactions Ar"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu metal");
    lmp->input->one("lattice fcc 4.08");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    lmp->input->one("kim_interactions Au");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // ASSERT_EQ(lmp->output->var_kim_periodic, 1);
    // TEST_FAILURE(".*ERROR: Incompatible units for KIM Simulator Model.*",
    //              lmp->input->one("kim_interactions Au"););


    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init LennardJones_Ar real");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    lmp->input->one("kim_interactions Ar");
    lmp->input->one("mass 1 39.95");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    int ifix = lmp->modify->find_fix("KIM_MODEL_STORE");
    ASSERT_GE(ifix, 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init LennardJones_Ar real");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    lmp->input->one("kim_interactions Ar");
    lmp->input->one("mass 1 39.95");
    lmp->input->one("run 1");
    lmp->input->one("kim_interactions Ar");
    lmp->input->one("run 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
}

TEST_F(KimCommandsTest, kim_param)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal kim_param command.*", lmp->input->one("kim_param"););
    TEST_FAILURE(".*ERROR: Incorrect arguments in kim_param command.\n"
                 "'kim_param get/set' is mandatory.*",
                 lmp->input->one("kim_param unknown shift 1 shift"););
    TEST_FAILURE(".*ERROR: Must use 'kim_init' before 'kim_param'.*",
                 lmp->input->one("kim_param get shift 1 shift"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu metal");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: kim_param can only be used with a KIM Portable Model.*",
                 lmp->input->one("kim_param get shift 1 shift"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init LennardJones612_UniversalShifted__MO_959249795837_003 real");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Illegal index '0' for "
                 "'shift' parameter with the extent of '1'.*",
                 lmp->input->one("kim_param get shift 0 shift"););
    TEST_FAILURE(".*ERROR: Illegal index '2' for "
                 "'shift' parameter with the extent of '1'.*",
                 lmp->input->one("kim_param get shift 2 shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range.\nExpected integer "
                 "parameter\\(s\\) instead of '1.' in index_range.*",
                 lmp->input->one("kim_param get shift 1. shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range '1-2' for 'shift' "
                 "parameter with the extent of '1'.*",
                 lmp->input->one("kim_param get shift 1:2 shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range.\nExpected integer "
                 "parameter\\(s\\) instead of '1-2' in index_range.*",
                 lmp->input->one("kim_param get shift 1-2 shift"););
    TEST_FAILURE(".*ERROR: Wrong number of arguments in 'kim_param "
                 "get' command.\nThe LAMMPS '3' variable names or "
                 "'s1 split' is mandatory.*",
                 lmp->input->one("kim_param get sigmas 1:3 s1 s2"););
    TEST_FAILURE(".*ERROR: Wrong argument in kim_param get command.\nThis "
                 "Model does not have the requested 'unknown' parameter.*",
                 lmp->input->one("kim_param get unknown 1 unknown"););
    TEST_FAILURE(".*ERROR: Wrong 'kim_param set' command.\n"
                 "To set the new parameter values, pair style must "
                 "be assigned.\nMust use 'kim_interactions' or"
                 "'pair_style kim' before 'kim_param set'.*",
                 lmp->input->one("kim_param set shift 1 2"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("kim_param get shift 1 shift");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_FALSE(lmp->input->variable->find("shift") == -1);
    ASSERT_TRUE(std::string(lmp->input->variable->retrieve("shift")) == std::string("1"));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init LennardJones612_UniversalShifted__MO_959249795837_003 real");
    lmp->input->one("lattice fcc 4.4300");
    lmp->input->one("region box block 0 10 0 10 0 10");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    lmp->input->one("kim_interactions Ar");
    lmp->input->one("mass 1 39.95");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Illegal index '2' for "
                 "'shift' parameter with the extent of '1'.*",
                 lmp->input->one("kim_param set shift 2 2"););
    TEST_FAILURE(".*ERROR: Illegal index_range.\nExpected integer "
                 "parameter\\(s\\) instead of '1.' in index_range.*",
                 lmp->input->one("kim_param set shift 1. shift"););
    TEST_FAILURE(".*ERROR: Illegal index_range '1-2' for "
                 "'shift' parameter with the extent of '1'.*",
                 lmp->input->one("kim_param set shift 1:2 2"););
    TEST_FAILURE(".*ERROR: Wrong number of variable values for pair coefficients.*",
                 lmp->input->one("kim_param set sigmas 1:3 0.5523570 0.4989030"););
    TEST_FAILURE(".*ERROR: Wrong argument for pair coefficients.\nThis "
                 "Model does not have the requested '0.4989030' parameter.*",
                 lmp->input->one("kim_param set sigmas 1:1 0.5523570 0.4989030"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("variable new_shift equal 2");
    lmp->input->one("kim_param set shift 1 ${new_shift}");
    lmp->input->one("kim_param get shift 1 shift");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_TRUE(std::string(lmp->input->variable->retrieve("shift")) == std::string("2"));
}

TEST_F(KimCommandsTest, kim_property)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();
    if (!LAMMPS::is_installed_pkg("PYTHON")) GTEST_SKIP();

    if (!lmp->python->has_minimum_version(3, 6)) {
        TEST_FAILURE(".*ERROR: Invalid Python version.\n"
                     "The kim-property Python package requires Python "
                     "3 >= 3.6 support.*",
                     lmp->input->one("kim_property"););
    } else {
        TEST_FAILURE(".*ERROR: Invalid kim_property command.*", 
                     lmp->input->one("kim_property"););
        TEST_FAILURE(".*ERROR: Invalid kim_property command.*",
                     lmp->input->one("kim_property create"););
        TEST_FAILURE(".*ERROR: Incorrect arguments in kim_property command.\n"
                     "'kim_property create/destroy/modify/remove/dump' "
                     "is mandatory.*",
                     lmp->input->one("kim_property unknown 1 atomic-mass"););
    }
#if defined(KIM_EXTRA_UNITTESTS)
        TEST_FAILURE(".*ERROR: Invalid 'kim_property create' command.*",
                     lmp->input->one("kim_property create 1"););
        TEST_FAILURE(".*ERROR: Invalid 'kim_property destroy' command.*",
                     lmp->input->one("kim_property destroy 1 cohesive-potential-energy-cubic-crystal"););
        TEST_FAILURE(".*ERROR: Invalid 'kim_property modify' command.*",
                     lmp->input->one("kim_property modify 1 key short-name"););
        TEST_FAILURE(".*ERROR: There is no property instance to modify the content.*",
                     lmp->input->one("kim_property modify 1 key short-name source-value 1 fcc"););
        TEST_FAILURE(".*ERROR: Invalid 'kim_property remove' command.*",
                     lmp->input->one("kim_property remove 1 key"););
        TEST_FAILURE(".*ERROR: There is no property instance to remove the content.*",
                     lmp->input->one("kim_property remove 1 key short-name"););
        TEST_FAILURE(".*ERROR: There is no property instance to dump the content.*",
                     lmp->input->one("kim_property dump results.edn"););
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp->input->one("clear");
        lmp->input->one("kim_init LennardJones612_UniversalShifted__MO_959249795837_003 real");
        lmp->input->one("kim_property create 1 cohesive-potential-energy-cubic-crystal");
        lmp->input->one("kim_property modify 1 key short-name source-value 1 fcc");
        lmp->input->one("kim_property destroy 1");
        if (!verbose) ::testing::internal::GetCapturedStdout();
#endif
}

TEST_F(KimCommandsTest, kim_query)
{
    if (!LAMMPS::is_installed_pkg("KIM")) GTEST_SKIP();
    
    TEST_FAILURE(".*ERROR: Illegal kim_query command.*", 
                 lmp->input->one("kim_query"););
    TEST_FAILURE(".*ERROR: Must use 'kim_init' before 'kim_query'.*",
                 lmp->input->one("kim_query a0 get_lattice_constant_cubic"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init LennardJones612_UniversalShifted__MO_959249795837_003 real");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Illegal kim_query command.\nThe keyword 'split' "
                 "must be followed by the name of the query function.*", 
                 lmp->input->one("kim_query a0 split"););

    TEST_FAILURE(".*ERROR: Illegal kim_query command.\nThe 'list' keyword "
                 "can not be used after 'split'.*", 
                 lmp->input->one("kim_query a0 split list"););

    TEST_FAILURE(".*ERROR: Illegal kim_query command.\nThe 'list' keyword "
                 "must be followed by \\('split' and\\) the name of the query "
                 "function.*", lmp->input->one("kim_query a0 list"););

    TEST_FAILURE(".*ERROR: Illegal 'model' key in kim_query command.*", 
                 lmp->input->one("kim_query a0 get_lattice_constant_cubic "
                                 "model=[MO_959249795837_003]"););

    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `crystal` "
                 "to kim_query is wrong. The query format is the "
                 "keyword=\\[value\\], where value is always an array of one "
                 "or more comma-separated items.*", 
                 lmp->input->one("kim_query a0 get_lattice_constant_cubic "
                                 "crystal"););

    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `"
                 "crystal=fcc` to kim_query is wrong. The query format is the "
                 "keyword=\\[value\\], where value is always an array of one "
                 "or more comma-separated items.*", 
                 lmp->input->one("kim_query a0 get_lattice_constant_cubic "
                                 "crystal=fcc"););

    TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `"
                 "crystal=\\[fcc` to kim_query is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*", 
                 lmp->input->one("kim_query a0 get_lattice_constant_cubic "
                                 "crystal=[fcc"););

   TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `"
                 "crystal=fcc\\]` to kim_query is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*", 
                 lmp->input->one("kim_query a0 get_lattice_constant_cubic "
                                 "crystal=fcc]"););
    
    std::string squery("kim_query a0 get_lattice_constant_cubic ");
    squery += "crystal=[\"fcc\"] species=\"Al\",\"Ni\" units=[\"angstrom\"]";

   TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "\"Al\",\"Ni\"` to kim_query is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*", 
                 lmp->input->one(squery););

    squery = "kim_query a0 get_lattice_constant_cubic ";
    squery += "crystal=[\"fcc\"] species=\"Al\",\"Ni\", units=[\"angstrom\"]";

   TEST_FAILURE(".*ERROR: Illegal query format.\nInput argument of `species="
                 "\"Al\",\"Ni\",` to kim_query is wrong. The query format is "
                 "the keyword=\\[value\\], where value is always an array of "
                 "one or more comma-separated items.*", 
                 lmp->input->one(squery););

    squery = "kim_query a0 get_lattice_constant_cubic crystal=[fcc] "
             "species=[Al]";
    TEST_FAILURE(".*ERROR: OpenKIM query failed:.*", lmp->input->one(squery););

    squery = "kim_query a0 get_lattice_constant_cubic crystal=[fcc] "
             "units=[\"angstrom\"]";
    TEST_FAILURE(".*ERROR: OpenKIM query failed:.*", lmp->input->one(squery););

#if defined(KIM_EXTRA_UNITTESTS)
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init EAM_Dynamo_Mendelev_2007_Zr__MO_848899341753_000 metal");
    
    squery = "kim_query latconst split get_lattice_constant_hexagonal ";
    squery += "crystal=[\"hcp\"] species=[\"Zr\"] units=[\"angstrom\"]";
    lmp->input->one(squery);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_TRUE((std::string(lmp->input->variable->retrieve("latconst_1")) == 
                 std::string("3.234055244384789")));
    ASSERT_TRUE((std::string(lmp->input->variable->retrieve("latconst_2")) == 
                 std::string("5.167650199630013")));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init EAM_Dynamo_Mendelev_2007_Zr__MO_848899341753_000 metal");
    
    squery = "kim_query latconst list get_lattice_constant_hexagonal ";
    squery += "crystal=[hcp] species=[Zr] units=[angstrom]";
    lmp->input->one(squery);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_TRUE((std::string(lmp->input->variable->retrieve("latconst")) == 
                 std::string("3.234055244384789  5.167650199630013")));

    squery = "kim_query latconst list get_lattice_constant_hexagonal ";
    squery += "crystal=[bcc] species=[Zr] units=[angstrom]";
    TEST_FAILURE(".*ERROR: OpenKIM query failed:.*", lmp->input->one(squery););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("kim_init EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal");
    
    squery = "kim_query alpha get_linear_thermal_expansion_coefficient_cubic ";
    squery += "crystal=[fcc] species=[Al] units=[1/K] temperature=[293.15] ";
    squery += "temperature_units=[K]";
    lmp->input->one(squery);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_TRUE((std::string(lmp->input->variable->retrieve("alpha")) == 
                 std::string("1.654960564704273e-05")));
#endif
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (have_openmpi && !LAMMPS_NS::Info::has_exceptions())
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
