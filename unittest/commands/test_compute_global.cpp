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

#include "../testing/core.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class ComputeGlobalTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ComputeGlobalTest";
        LAMMPSTest::SetUp();
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in.fourmol\"");
            command("group allwater molecule 3:6");
            command("region half block 0.0 INF INF INF INF INF");
            END_HIDE_OUTPUT();
        }
    }

    double get_scalar(const char *id)
    {
        return *(double *)lammps_extract_compute(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);
    }

    double *get_vector(const char *id)
    {
        return (double *)lammps_extract_compute(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    }

    double **get_array(const char *id)
    {
        return (double **)lammps_extract_compute(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY);
    }
};

TEST_F(ComputeGlobalTest, Energy)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();
    int has_tally = lammps_config_has_package("TALLY");

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut/coul/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");

    command("compute ke1 all ke");
    command("compute ke2 allwater ke");
    command("compute pe1 all pe");
    command("compute pe2 all pe bond");
    command("compute pe3 all pe angle dihedral");
    command("compute pr1 all pressure thermo_temp");
    command("compute pr2 all pressure NULL virial");
    command("compute pr3 all pressure NULL angle dihedral");
    std::string thermo_style = "c_ke1 c_ke2 c_pe1 c_pe2 c_pe3 c_pr1 c_pr2 c_pr3";

    if (has_tally) {
        command("compute pe4 all pe/tally allwater");
        command("compute pe5 all pe/mol/tally all");
        command("compute pe6 all pe pair");
        thermo_style += " c_pe4 c_pe5[*]";
    }

    command("thermo_style custom " + thermo_style);
    command("run 0 post no");
    END_HIDE_OUTPUT();

    EXPECT_DOUBLE_EQ(get_scalar("ke1"), 2.3405256449146168);
    EXPECT_DOUBLE_EQ(get_scalar("ke2"), 1.192924237073665);
    EXPECT_DOUBLE_EQ(get_scalar("pe1"), 24155.155261642241);
    EXPECT_DOUBLE_EQ(get_scalar("pe2"), 361.37528652881286);
    EXPECT_DOUBLE_EQ(get_scalar("pe3"), 0.0);
    EXPECT_NEAR(get_scalar("pr1"), 1956948.4735454607, 0.000000005);
    EXPECT_NEAR(get_scalar("pr2"), 1956916.7725807722, 0.000000005);
    EXPECT_DOUBLE_EQ(get_scalar("pr3"), 0.0);
    auto pr1 = get_vector("pr1");
    auto pr2 = get_vector("pr2");
    auto pr3 = get_vector("pr3");
    EXPECT_NEAR(pr1[0], 2150600.9207200543, 0.000000005);
    EXPECT_NEAR(pr1[1], 1466949.7512112649, 0.000000005);
    EXPECT_NEAR(pr1[2], 2253294.7487050635, 0.000000005);
    EXPECT_NEAR(pr1[3], 856643.16926486336, 0.000000005);
    EXPECT_NEAR(pr1[4], 692710.86929464422, 0.000000005);
    EXPECT_NEAR(pr1[5], -44403.909298603547, 0.000000005);
    EXPECT_NEAR(pr2[0], 2150575.6989334146, 0.000000005);
    EXPECT_NEAR(pr2[1], 1466911.3911461537, 0.000000005);
    EXPECT_NEAR(pr2[2], 2253263.2276627473, 0.000000005);
    EXPECT_NEAR(pr2[3], 856632.34707690508, 0.000000005);
    EXPECT_NEAR(pr2[4], 692712.89222328411, 0.000000005);
    EXPECT_NEAR(pr2[5], -44399.277068014424, 0.000000005);
    EXPECT_DOUBLE_EQ(pr3[0], 0.0);
    EXPECT_DOUBLE_EQ(pr3[1], 0.0);
    EXPECT_DOUBLE_EQ(pr3[2], 0.0);
    EXPECT_DOUBLE_EQ(pr3[3], 0.0);
    EXPECT_DOUBLE_EQ(pr3[4], 0.0);
    EXPECT_DOUBLE_EQ(pr3[5], 0.0);

    if (has_tally) {
        EXPECT_NEAR(get_scalar("pe4"), 15425.840923850392, 0.000000005);
        auto pe5 = get_vector("pe5");
        EXPECT_NEAR(pe5[0], 23803.966677151559, 0.000000005);
        EXPECT_NEAR(pe5[1], -94.210004432380643, 0.000000005);
        EXPECT_NEAR(pe5[2], 115.58040355478101, 0.000000005);
        EXPECT_NEAR(pe5[3], -31.557101160514257, 0.000000005);
    }

    TEST_FAILURE(".*ERROR: Compute pressure must use group all.*",
                 command("compute pr5 allwater pressure thermo_temp"););
    TEST_FAILURE(".*ERROR: Compute pressure requires temperature ID to include kinetic energy.*",
                 command("compute pr5 all pressure NULL"););
    TEST_FAILURE(".*ERROR: Could not find compute pressure temperature ID",
                 command("compute pr5 all pressure xxx"););

    TEST_FAILURE(".*ERROR: Reuse of compute ID 'pe2'.*", command("compute pe2 all pe"););
    TEST_FAILURE(".*ERROR: Compute pe must use group all.*", command("compute pe allwater pe"););
    TEST_FAILURE(".*ERROR: Illegal compute command.*", command("compute pe potential"););
}

TEST_F(ComputeGlobalTest, Geometry)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();
    int has_extra = lammps_config_has_package("EXTRA-COMPUTE");

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");

    command("compute com1 all com");
    command("compute com2 allwater com");
    command("compute mu1 all dipole");
    command("compute mu2 allwater dipole geometry ");
    command("compute rg1 all gyration");
    command("compute rg2 allwater gyration");
    std::string thermo_style = "c_com1[*] c_com2[*] c_rg1[*] c_rg2[*]";

    if (has_extra) {
        command("compute mom1 all momentum");
        command("compute mom2 allwater momentum");
        command("compute mop1 all stress/mop x 0.0 total");
        command("compute mop2 all stress/mop/profile z lower 0.5 kin conf");
        thermo_style += " c_mu1 c_mu2 c_mop1[*] c_mop2[1][1]";
    }

    command("thermo_style custom " + thermo_style);
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto com1 = get_vector("com1");
    auto com2 = get_vector("com2");
    auto mu1  = get_vector("mu1");
    auto mu2  = get_vector("mu2");
    auto rg1  = get_vector("rg1");
    auto rg2  = get_vector("rg2");

    EXPECT_NEAR(com1[0], 1.4300952724948282, 0.0000000005);
    EXPECT_NEAR(com1[1], -0.29759806705328351, 0.0000000005);
    EXPECT_NEAR(com1[2], -0.7245120195899285, 0.0000000005);
    EXPECT_NEAR(com2[0], 1.7850913321989679, 0.0000000005);
    EXPECT_NEAR(com2[1], -0.45168408952146238, 0.0000000005);
    EXPECT_NEAR(com2[2], -0.60215022088294912, 0.0000000005);

    EXPECT_NEAR(get_scalar("mu1"), 1.8335537504770163, 0.0000000005);
    EXPECT_NEAR(get_scalar("mu2"), 1.7849382239204072, 0.0000000005);
    EXPECT_NEAR(mu1[0], 0.41613191281297729, 0.0000000005);
    EXPECT_NEAR(mu1[1], 1.0056523085627747, 0.0000000005);
    EXPECT_NEAR(mu1[2], -1.4756073398127658, 0.0000000005);
    EXPECT_NEAR(mu2[0], -0.029474795088977768, 0.0000000005);
    EXPECT_NEAR(mu2[1], 1.153516133030746, 0.0000000005);
    EXPECT_NEAR(mu2[2], -1.3618135814069394, 0.0000000005);

    EXPECT_NEAR(get_scalar("rg1"), 3.8495643473797196, 0.0000000005);
    EXPECT_NEAR(get_scalar("rg2"), 5.4558163385611342, 0.0000000005);
    EXPECT_NEAR(rg1[0], 3.6747807397432752, 0.0000000005);
    EXPECT_NEAR(rg1[1], 6.5440303159316278, 0.0000000005);
    EXPECT_NEAR(rg1[2], 4.6003346089421457, 0.0000000005);
    EXPECT_NEAR(rg1[3], -0.4639249501367636, 0.0000000005);
    EXPECT_NEAR(rg1[4], -1.8859032304357459, 0.0000000005);
    EXPECT_NEAR(rg1[5], 0.2339161878440186, 0.0000000005);
    EXPECT_NEAR(rg2[0], 6.2582260148310143, 0.0000000005);
    EXPECT_NEAR(rg2[1], 13.353763805454184, 0.0000000005);
    EXPECT_NEAR(rg2[2], 10.153942099825425, 0.0000000005);
    EXPECT_NEAR(rg2[3], 1.2965604701522486, 0.0000000005);
    EXPECT_NEAR(rg2[4], -5.0315240817290841, 0.0000000005);
    EXPECT_NEAR(rg2[5], 1.1103378503822141, 0.0000000005);
    if (has_extra) {
        auto mom1 = get_vector("mom1");
        auto mom2 = get_vector("mom2");
        auto mop1 = get_vector("mop1");
        auto mop2 = get_array("mop2");
        EXPECT_DOUBLE_EQ(mom1[0], 0.0054219056685341164);
        EXPECT_DOUBLE_EQ(mom1[1], -0.054897225112275558);
        EXPECT_DOUBLE_EQ(mom1[2], 0.059097392692385661);
        EXPECT_DOUBLE_EQ(mom2[0], -0.022332069630161717);
        EXPECT_DOUBLE_EQ(mom2[1], -0.056896553865696115);
        EXPECT_DOUBLE_EQ(mom2[2], 0.069179891052881484);
        EXPECT_DOUBLE_EQ(mop1[0], 3522311.3572200728);
        EXPECT_DOUBLE_EQ(mop1[1], 2871104.9055934539);
        EXPECT_DOUBLE_EQ(mop1[2], -4136077.5224247416);
        EXPECT_DOUBLE_EQ(mop2[0][0], -8.0869239999999998);
        EXPECT_DOUBLE_EQ(mop2[0][1], 0.0);
        EXPECT_DOUBLE_EQ(mop2[0][2], 0.0);
        EXPECT_DOUBLE_EQ(mop2[1][0], -7.5869239999999998);
        EXPECT_DOUBLE_EQ(mop2[1][1], 0.0);
        EXPECT_DOUBLE_EQ(mop2[1][2], 0.0);
    }
}

TEST_F(ComputeGlobalTest, Reduction)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");

    command("variable v atom sqrt(vx*vx+vy*vy+vz*vz)");
    command("variable id atom id");
    command("fix chg all store/state 0 q");
    command("compute ke all ke/atom");
    command("compute min allwater reduce min x fx v_v");
    command("compute chg all reduce max f_chg");
    command("compute max all reduce max y fy v_v");
    command("compute ave all reduce/region half ave z fz v_v");
    command("compute sum allwater reduce/region half sum vx vy vz");
    command("compute rep all reduce max v_id v_v v_id y replace 1 2 replace 3 4");
    std::string thermo_style = "c_min[*] c_chg c_max[*] c_sum[*] c_ave[*] c_rep[*]";

    command("thermo_style custom " + thermo_style);
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto min = get_vector("min");
    auto max = get_vector("max");
    auto sum = get_vector("sum");
    auto ave = get_vector("ave");
    auto rep = get_vector("rep");

    EXPECT_DOUBLE_EQ(get_scalar("chg"), 0.51000000000000001);

    EXPECT_DOUBLE_EQ(min[0], -2.7406520384725965);
    EXPECT_DOUBLE_EQ(min[1], -20385.448391361348);
    EXPECT_DOUBLE_EQ(min[2], 0.00071995632406981081);

    EXPECT_DOUBLE_EQ(max[0], 4.0120175892854135);
    EXPECT_DOUBLE_EQ(max[1], 21193.39005673242);
    EXPECT_DOUBLE_EQ(max[2], 0.0072167889062371513);

    EXPECT_DOUBLE_EQ(sum[0], 0.0021436162503408024);
    EXPECT_DOUBLE_EQ(sum[1], -0.013760203913131267);
    EXPECT_DOUBLE_EQ(sum[2], 0.017517003988402391);

    EXPECT_DOUBLE_EQ(ave[0], -1.3013763067943667);
    EXPECT_DOUBLE_EQ(ave[1], -619.60864441905312);
    EXPECT_DOUBLE_EQ(ave[2], 0.0035263629500884397);

    // index of max v_v
    EXPECT_DOUBLE_EQ(rep[0], 20);
    EXPECT_DOUBLE_EQ(rep[1], max[2]);
    // index of max y
    EXPECT_DOUBLE_EQ(rep[2], 26);
    EXPECT_DOUBLE_EQ(rep[3], max[0]);
}

TEST_F(ComputeGlobalTest, Counts)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style zero 10.0");
    command("pair_coeff * *");

    command("variable t1 atom type==1");
    command("variable t2 atom type==2");
    command("variable t3 atom type==3");
    command("variable t4 atom type==4");
    command("variable t5 atom type==5");
    command("compute tsum all reduce sum v_t1 v_t2 v_t3 v_t4 v_t5");
    command("compute tcnt all count/type atom");
    command("compute bcnt all count/type bond");
    command("compute acnt all count/type angle");
    command("compute dcnt all count/type dihedral");
    command("compute icnt all count/type improper");
    command("thermo_style custom c_tsum[*] c_tcnt[*] c_bcnt[*] c_acnt[*] c_dcnt[*] c_icnt[*]");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto tsum = get_vector("tsum");
    auto tcnt = get_vector("tcnt");
    auto bcnt = get_vector("bcnt");
    auto bbrk = get_scalar("bcnt");
    auto acnt = get_vector("acnt");
    auto dcnt = get_vector("dcnt");
    auto icnt = get_vector("icnt");

    EXPECT_DOUBLE_EQ(tsum[0], tcnt[0]);
    EXPECT_DOUBLE_EQ(tsum[1], tcnt[1]);
    EXPECT_DOUBLE_EQ(tsum[2], tcnt[2]);
    EXPECT_DOUBLE_EQ(tsum[3], tcnt[3]);
    EXPECT_DOUBLE_EQ(tsum[4], tcnt[4]);

    EXPECT_DOUBLE_EQ(bbrk, 0.0);

    EXPECT_DOUBLE_EQ(bcnt[0], 3.0);
    EXPECT_DOUBLE_EQ(bcnt[1], 6.0);
    EXPECT_DOUBLE_EQ(bcnt[2], 3.0);
    EXPECT_DOUBLE_EQ(bcnt[3], 2.0);
    EXPECT_DOUBLE_EQ(bcnt[4], 10.0);

    EXPECT_DOUBLE_EQ(acnt[0], 6.0);
    EXPECT_DOUBLE_EQ(acnt[1], 10.0);
    EXPECT_DOUBLE_EQ(acnt[2], 5.0);
    EXPECT_DOUBLE_EQ(acnt[3], 9.0);

    EXPECT_DOUBLE_EQ(dcnt[0], 3.0);
    EXPECT_DOUBLE_EQ(dcnt[1], 8.0);
    EXPECT_DOUBLE_EQ(dcnt[2], 3.0);
    EXPECT_DOUBLE_EQ(dcnt[3], 4.0);
    EXPECT_DOUBLE_EQ(dcnt[4], 13.0);

    EXPECT_DOUBLE_EQ(icnt[0], 1.0);
    EXPECT_DOUBLE_EQ(icnt[1], 1.0);

    BEGIN_HIDE_OUTPUT();
    command("delete_bonds all bond 3 remove");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    bcnt = get_vector("bcnt");
    bbrk = get_scalar("bcnt");
    acnt = get_vector("acnt");
    dcnt = get_vector("dcnt");
    icnt = get_vector("icnt");

    EXPECT_DOUBLE_EQ(bbrk, 0.0);
    EXPECT_DOUBLE_EQ(bcnt[0], 3.0);
    EXPECT_DOUBLE_EQ(bcnt[1], 6.0);
    EXPECT_DOUBLE_EQ(bcnt[2], 0.0);
    EXPECT_DOUBLE_EQ(bcnt[3], 2.0);
    EXPECT_DOUBLE_EQ(bcnt[4], 10.0);

    EXPECT_DOUBLE_EQ(acnt[0], 6.0);
    EXPECT_DOUBLE_EQ(acnt[1], 10.0);
    EXPECT_DOUBLE_EQ(acnt[2], 5.0);
    EXPECT_DOUBLE_EQ(acnt[3], 9.0);

    EXPECT_DOUBLE_EQ(dcnt[0], 3.0);
    EXPECT_DOUBLE_EQ(dcnt[1], 8.0);
    EXPECT_DOUBLE_EQ(dcnt[2], 3.0);
    EXPECT_DOUBLE_EQ(dcnt[3], 4.0);
    EXPECT_DOUBLE_EQ(dcnt[4], 13.0);

    EXPECT_DOUBLE_EQ(icnt[0], 1.0);
    EXPECT_DOUBLE_EQ(icnt[1], 1.0);
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (LAMMPS_NS::platform::mpi_vendor() == "Open MPI" && !Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. Death tests will be skipped\n";

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
