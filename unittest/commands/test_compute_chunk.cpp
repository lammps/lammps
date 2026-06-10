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
#include "lammps.h"
#include "library.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

// we compare floating point numbers with 8 digits precision after the decimal point
static constexpr double EPSILON = 1.0e-8;

namespace LAMMPS_NS {

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class ComputeChunkTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ComputeChunkTest";
        LAMMPSTest::SetUp();
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in.fourmol\"");
            command("group allwater molecule 3:6");
            command("region half block 0.0 INF INF INF INF INF");
            command("compute tags all property/atom id");
            command("compute bin1d all chunk/atom bin/1d x lower 3.0 units box");
            command("compute bin2d all chunk/atom bin/2d x lower 3.0 y lower 3.0 units box");
            command("compute bin3d all chunk/atom bin/3d x lower 3.0 y lower 3.0 z lower 3.0 units "
                    "box");
            command("compute binsph all chunk/atom bin/sphere 0.0 0.0 0.0 0.01 6.01 6 units box");
            command("compute bincyl all chunk/atom bin/cylinder z lower 3.0 1.0 1.0 0.01 6.01 6 "
                    "units box");
            command("compute mols all chunk/atom molecule");
            command("compute types all chunk/atom type");
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

    double *get_peratom(const char *id)
    {
        return (double *)lammps_extract_compute(lmp, id, LMP_STYLE_ATOM, LMP_TYPE_VECTOR);
    }
};

static constexpr int chunk1d[]  = {0, 2, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 3,
                                   4, 3, 3, 3, 3, 3, 4, 4, 4, 3, 3, 3, 2, 2, 2};
static constexpr int chalf1d[]  = {0, 0, 3, 0, 0, 0, 3, 3, 3, 3, 3, 3, 4, 4, 3,
                                   4, 3, 3, 3, 3, 3, 4, 4, 4, 3, 3, 3, 0, 0, 0};
static constexpr int chunk2d[]  = {0,  9,  14, 8,  9,  8,  13, 13, 13, 13, 13, 12, 18, 18, 13,
                                   18, 12, 12, 14, 14, 14, 17, 17, 17, 14, 14, 14, 7,  7,  7};
static constexpr int chunk3d[]  = {0,  43, 68, 38, 43, 38, 63, 62, 63, 63, 63, 58, 88, 88, 62,
                                   88, 58, 59, 67, 67, 67, 82, 82, 82, 69, 69, 69, 34, 34, 34};
static constexpr int chunksph[] = {0, 3, 4, 2, 3, 2, 2, 3, 2, 2, 3, 4, 4, 5, 4,
                                   4, 4, 4, 6, 6, 6, 6, 6, 6, 5, 5, 6, 6, 6, 5};
static constexpr int chunkcyl[] = {0,  8,  13, 8,  13, 8,  8,  7,  8,  8,  13, 18, 13, 18, 12,
                                   13, 18, 19, 12, 7,  17, 27, 27, 27, 14, 14, 19, 29, 29, 29};
static constexpr int chunkmol[] = {0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
                                   2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6};
static constexpr int chunktyp[] = {0, 3, 2, 1, 2, 2, 1, 4, 3, 2, 1, 2, 1, 2, 2,
                                   2, 1, 4, 4, 2, 2, 5, 2, 2, 5, 2, 2, 5, 2, 2};

TEST_F(ComputeChunkTest, ChunkAtom)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut/coul/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");
    command("dump 1 all custom 1 compute_chunk_atom.lammpstrj "
            "id c_bin1d c_bin2d c_bin3d c_binsph c_bincyl c_mols c_types c_tags");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int natoms = lammps_get_natoms(lmp);
    EXPECT_EQ(get_scalar("bin1d"), 5);
    EXPECT_EQ(get_scalar("bin2d"), 25);
    EXPECT_EQ(get_scalar("bin3d"), 125);
    EXPECT_EQ(get_scalar("binsph"), 6);
    EXPECT_EQ(get_scalar("bincyl"), 30);
    EXPECT_EQ(get_scalar("mols"), 6);
    EXPECT_EQ(get_scalar("types"), 5);

    auto *cbin1d  = get_peratom("bin1d");
    auto *cbin2d  = get_peratom("bin2d");
    auto *cbin3d  = get_peratom("bin3d");
    auto *cbinsph = get_peratom("binsph");
    auto *cbincyl = get_peratom("bincyl");
    auto *cmols   = get_peratom("mols");
    auto *ctypes  = get_peratom("types");
    auto *tag     = get_peratom("tags");

    for (int i = 0; i < natoms; ++i) {
        EXPECT_EQ(cbin1d[i], chunk1d[(int)tag[i]]);
        EXPECT_EQ(cbin2d[i], chunk2d[(int)tag[i]]);
        EXPECT_EQ(cbin3d[i], chunk3d[(int)tag[i]]);
        EXPECT_EQ(cbinsph[i], chunksph[(int)tag[i]]);
        EXPECT_EQ(cbincyl[i], chunkcyl[(int)tag[i]]);
        EXPECT_EQ(cmols[i], chunkmol[(int)tag[i]]);
        EXPECT_EQ(ctypes[i], chunktyp[(int)tag[i]]);
    }

    BEGIN_HIDE_OUTPUT();
    command("uncompute bin1d");
    command("compute bin1d all chunk/atom bin/1d x lower 0.2 units reduced region half");
    command("uncompute bin3d");
    command("compute bin3d all chunk/atom bin/3d x lower 3.0 y lower 3.0 z lower 3.0 "
            "compress yes units box");
    END_HIDE_OUTPUT();
    EXPECT_EQ(get_scalar("bin1d"), 5);
    EXPECT_EQ(get_scalar("bin3d"), 12);

    cbin1d = get_peratom("bin1d");
    tag    = get_peratom("tags");
    for (int i = 0; i < natoms; ++i) {
        EXPECT_EQ(cbin1d[i], chalf1d[(int)tag[i]]);
    }

    // cleanup
    platform::unlink("compute_chunk_atom.lammpstrj");
}

TEST_F(ComputeChunkTest, PropertyChunk)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut/coul/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");
    command("uncompute bin3d");
    command("compute bin3d all chunk/atom bin/3d x lower 3.0 y lower 3.0 z lower 3.0 "
            "compress yes units box");
    command("compute prop1 all property/chunk bin1d count");
    command("compute prop2 all property/chunk bin2d count");
    command("compute prop3 all property/chunk bin3d id count");
    command("fix hist1 all ave/time 1 1 1 c_prop1 mode vector");
    command("fix hist2 all ave/time 1 1 1 c_prop2 mode vector");
    command("fix hist3 all ave/time 1 1 1 c_prop3[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *cprop1 = get_vector("prop1");
    EXPECT_EQ(cprop1[0], 0);
    EXPECT_EQ(cprop1[1], 7);
    EXPECT_EQ(cprop1[2], 16);
    EXPECT_EQ(cprop1[3], 6);
    EXPECT_EQ(cprop1[4], 0);

    auto *cprop2 = get_vector("prop2");
    int nempty   = 0;
    int ncount   = 0;
    for (int i = 0; i < 25; ++i) {
        if (cprop2[i] == 0)
            ++nempty;
        else
            ncount += cprop2[i];
    }
    EXPECT_EQ(nempty, 17);
    EXPECT_EQ(ncount, 29);

    auto *cprop3 = get_array("prop3");
    EXPECT_EQ(cprop3[0][0], 34);
    EXPECT_EQ(cprop3[1][0], 38);
    EXPECT_EQ(cprop3[2][0], 43);
    EXPECT_EQ(cprop3[3][0], 58);
    EXPECT_EQ(cprop3[4][0], 59);
    EXPECT_EQ(cprop3[5][0], 62);
    EXPECT_EQ(cprop3[6][0], 63);
    EXPECT_EQ(cprop3[7][0], 67);
    EXPECT_EQ(cprop3[8][0], 68);
    EXPECT_EQ(cprop3[9][0], 69);
    EXPECT_EQ(cprop3[10][0], 82);
    EXPECT_EQ(cprop3[11][0], 88);

    EXPECT_EQ(cprop3[0][1], 3);
    EXPECT_EQ(cprop3[1][1], 2);
    EXPECT_EQ(cprop3[2][1], 2);
    EXPECT_EQ(cprop3[3][1], 2);
    EXPECT_EQ(cprop3[4][1], 1);
    EXPECT_EQ(cprop3[5][1], 2);
    EXPECT_EQ(cprop3[6][1], 4);
    EXPECT_EQ(cprop3[7][1], 3);
    EXPECT_EQ(cprop3[8][1], 1);
    EXPECT_EQ(cprop3[9][1], 3);
    EXPECT_EQ(cprop3[10][1], 3);
    EXPECT_EQ(cprop3[11][1], 3);
}

TEST_F(ComputeChunkTest, ChunkComputes)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut/coul/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");
    command("compute ang all angmom/chunk mols");
    command("compute com all com/chunk mols");
    command("compute dip all dipole/chunk mols geometry");
    command("compute gyr all gyration/chunk mols");
    command("compute mom all inertia/chunk mols");
    command("compute omg all omega/chunk mols");
    command("compute tmp all temp/chunk mols com yes");
    command("compute trq all torque/chunk mols");
    command("compute vcm all vcm/chunk mols");
    command("fix hist1 all ave/time 1 1 1 c_ang[*] c_com[*] c_dip[*] c_gyr c_mom[*] c_omg[*] "
            "c_trq[*] c_vcm[*] mode vector");
    command("fix hist2 all ave/time 1 1 1 c_tmp mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    auto *cang = get_array("ang");
    auto *ccom = get_array("com");
    auto *cdip = get_array("dip");
    auto *cgyr = get_vector("gyr");
    auto *cmom = get_array("mom");
    auto *comg = get_array("omg");
    auto *ctmp = get_vector("tmp");
    auto *ctrq = get_array("trq");
    auto *cvcm = get_array("vcm");
    EXPECT_NEAR(cang[0][0], -0.01906982, EPSILON);
    EXPECT_NEAR(cang[0][1], -0.02814532, EPSILON);
    EXPECT_NEAR(cang[0][2], -0.03357393, EPSILON);
    EXPECT_NEAR(cang[5][0], 0.00767837, EPSILON);
    EXPECT_NEAR(cang[5][1], -0.00303138, EPSILON);
    EXPECT_NEAR(cang[5][2], 0.00740977, EPSILON);
    EXPECT_NEAR(ccom[1][0], 2.27051374, EPSILON);
    EXPECT_NEAR(ccom[1][1], -1.21038876, EPSILON);
    EXPECT_NEAR(ccom[1][2], -0.58581655, EPSILON);
    EXPECT_NEAR(ccom[5][0], -1.98284693, EPSILON);
    EXPECT_NEAR(ccom[5][1], -4.17351226, EPSILON);
    EXPECT_NEAR(ccom[5][2], 2.04850072, EPSILON);
    EXPECT_NEAR(cmom[2][0], 4.28810281, EPSILON);
    EXPECT_NEAR(cmom[2][1], 4.99562488, EPSILON);
    EXPECT_NEAR(cmom[2][2], 5.34954800, EPSILON);
    EXPECT_NEAR(cmom[5][0], 3.06867233, EPSILON);
    EXPECT_NEAR(cmom[5][1], 5.24202887, EPSILON);
    EXPECT_NEAR(cmom[5][2], 6.06478557, EPSILON);
    EXPECT_NEAR(comg[3][0], -0.00349423, EPSILON);
    EXPECT_NEAR(comg[3][1], -0.00025062, EPSILON);
    EXPECT_NEAR(comg[3][2], -0.00323573, EPSILON);
    EXPECT_NEAR(comg[5][0], 0.00437315, EPSILON);
    EXPECT_NEAR(comg[5][1], 0.00029335, EPSILON);
    EXPECT_NEAR(comg[5][2], 0.00268517, EPSILON);
    EXPECT_NEAR(ctrq[4][0], -0.94086086, EPSILON);
    EXPECT_NEAR(ctrq[4][1], 0.56227336, EPSILON);
    EXPECT_NEAR(ctrq[4][2], 0.75139995, EPSILON);
    EXPECT_NEAR(ctrq[5][0], -0.07066910, EPSILON);
    EXPECT_NEAR(ctrq[5][1], -0.58556032, EPSILON);
    EXPECT_NEAR(ctrq[5][2], -0.81513604, EPSILON);
    EXPECT_NEAR(cvcm[0][0], -0.00011274, EPSILON);
    EXPECT_NEAR(cvcm[0][1], 0.00071452, EPSILON);
    EXPECT_NEAR(cvcm[0][2], -0.00017908, EPSILON);
    EXPECT_NEAR(cvcm[5][0], -0.00063326, EPSILON);
    EXPECT_NEAR(cvcm[5][1], 0.00007092, EPSILON);
    EXPECT_NEAR(cvcm[5][2], 0.00045545, EPSILON);
    EXPECT_NEAR(cdip[0][3], 0.35912150, EPSILON);
    EXPECT_NEAR(cdip[1][3], 0.68453713, EPSILON);
    EXPECT_NEAR(cdip[2][3], 0.50272643, EPSILON);
    EXPECT_NEAR(cdip[3][3], 0.50845862, EPSILON);
    EXPECT_NEAR(cdip[4][3], 0.49757365, EPSILON);
    EXPECT_NEAR(cdip[5][3], 0.49105019, EPSILON);
    EXPECT_NEAR(cgyr[0], 1.48351858, EPSILON);
    EXPECT_NEAR(cgyr[1], 1.56649567, EPSILON);
    EXPECT_NEAR(cgyr[2], 0.55196552, EPSILON);
    EXPECT_NEAR(cgyr[3], 0.54573649, EPSILON);
    EXPECT_NEAR(cgyr[4], 0.54793875, EPSILON);
    EXPECT_NEAR(cgyr[5], 0.54708204, EPSILON);
    EXPECT_NEAR(ctmp[0], 1.08268576, EPSILON);
    EXPECT_NEAR(ctmp[1], 1.61905718, EPSILON);
    EXPECT_NEAR(ctmp[2], 1.41991778, EPSILON);
    EXPECT_NEAR(ctmp[3], 0.55484671, EPSILON);
    EXPECT_NEAR(ctmp[4], -0.06062938, EPSILON);
    EXPECT_NEAR(ctmp[5], -0.09219489, EPSILON);
}

TEST_F(ComputeChunkTest, ChunkTIP4PComputes)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();
    if (!info->has_style("compute", "dipole/tip4p/chunk")) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style tip4p/cut 5 2 5 1 0.15 10.0");
    command("pair_coeff * *");
    command("bond_coeff * 0.9572");
    command("angle_coeff * 104.52");
    command("compute dip all dipole/tip4p/chunk mols geometry");
    command("fix hist1 all ave/time 1 1 1 c_dip[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    auto *cdip = get_array("dip");
    EXPECT_NEAR(cdip[0][3], 0.35912150, EPSILON);
    EXPECT_NEAR(cdip[1][3], 0.68453713, EPSILON);
    EXPECT_NEAR(cdip[2][3], 0.50272643, EPSILON);
    EXPECT_NEAR(cdip[3][3], 0.37828094, EPSILON);
    EXPECT_NEAR(cdip[4][3], 0.37018279, EPSILON);
    EXPECT_NEAR(cdip[5][3], 0.36532949, EPSILON);
}

TEST_F(ComputeChunkTest, ChunkSpreadGlobal)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut/coul/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");
    command("compute gyr all gyration/chunk mols");
    command("compute spr all chunk/spread/atom mols c_gyr");
    command("compute glb all global/atom c_mols c_gyr");
    command("variable odd atom ((id+1)%2)+1");
    command("compute odd all global/atom v_odd c_gyr");
    command("fix ave all ave/atom 1 1 1 c_spr c_glb c_odd");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int natoms = lammps_get_natoms(lmp);

    auto *cgyr = get_vector("gyr");
    auto *cspr = get_peratom("spr");
    auto *cglb = get_peratom("glb");
    auto *codd = get_peratom("odd");
    auto *ctag = get_peratom("tags");

    for (int i = 0; i < natoms; ++i) {
        EXPECT_EQ(cspr[i], cgyr[chunkmol[(int)ctag[i]] - 1]);
        EXPECT_EQ(cglb[i], cgyr[chunkmol[(int)ctag[i]] - 1]);
        EXPECT_EQ(codd[i], cgyr[(((int)ctag[i] + 1) % 2)]);
    }
}

TEST_F(ComputeChunkTest, ChunkReduce)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut/coul/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");
    command("compute prp all property/chunk mols count");
    command("variable one atom 1");
    command("compute red all reduce/chunk mols sum v_one");
    command("fix ave all ave/time 1 1 1 c_prp c_red mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int nchunks = get_scalar("mols");

    auto *cprp = get_vector("prp");
    auto *cred = get_vector("red");

    for (int i = 0; i < nchunks; ++i)
        EXPECT_EQ(cprp[i], cred[i]);
}

// finite-size particles must contribute their own moment of inertia
// to compute inertia/chunk (issue #3710).

class ComputeInertiaChunkTest : public ComputeChunkTest {
protected:
    void SetUp() override
    {
        testbinary = "ComputeInertiaChunkTest";
        LAMMPSTest::SetUp();
    }
};

TEST_F(ComputeInertiaChunkTest, Ellipsoid)
{
    if (!info->has_style("atom", "ellipsoid")) GTEST_SKIP();

    // single chunk = both atoms (type 1); same configuration and analytic
    // result as the ComputeInertiaTest.Ellipsoid global-compute test:
    // two axis-aligned ellipsoids, semi-axes a=1,b=2,c=3, m1=1 at origin,
    // m2=3 at (5,0,0) -> tensor (10.4, 26.75, 22.75, 0, 0, 0)

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style ellipsoid");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 5.0 0.0 0.0 units box");
    command("set group all shape 2.0 4.0 6.0");
    command("set atom 1 mass 1.0");
    command("set atom 2 mass 3.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("compute ch all chunk/atom type");
    command("compute icc all inertia/chunk ch");
    command("fix h all ave/time 1 1 1 c_icc[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *icc = get_array("icc");
    EXPECT_NEAR(icc[0][0], 10.4, 1.0e-12);
    EXPECT_NEAR(icc[0][1], 26.75, 1.0e-12);
    EXPECT_NEAR(icc[0][2], 22.75, 1.0e-12);
    EXPECT_NEAR(icc[0][3], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][4], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][5], 0.0, 1.0e-12);
}

TEST_F(ComputeInertiaChunkTest, Sphere)
{
    if (!info->has_style("atom", "sphere")) GTEST_SKIP();

    // single chunk = both atoms (type 1); two finite spheres of radius 1,
    // mass 1, at (0,0,0) and (4,0,0) -> tensor (0.8, 8.8, 8.8, 0, 0, 0)

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style sphere");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 4.0 0.0 0.0 units box");
    command("set group all diameter 2.0");
    command("set group all mass 1.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("compute ch all chunk/atom type");
    command("compute icc all inertia/chunk ch");
    command("fix h all ave/time 1 1 1 c_icc[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *icc = get_array("icc");
    EXPECT_NEAR(icc[0][0], 0.8, 1.0e-12);
    EXPECT_NEAR(icc[0][1], 8.8, 1.0e-12);
    EXPECT_NEAR(icc[0][2], 8.8, 1.0e-12);
    EXPECT_NEAR(icc[0][3], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][4], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][5], 0.0, 1.0e-12);
}

TEST_F(ComputeInertiaChunkTest, Superellipsoid)
{
    if (!info->has_style("atom", "ellipsoid")) GTEST_SKIP();

    // same as Ellipsoid, but with atom_style "ellipsoid superellipsoid".
    // With default blockiness (2,2) the result is identical and the
    // superellipsoid (bonus_super) branch is exercised.  Mass must be set
    // before the shape, since the stored inertia is computed at set-shape.

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style ellipsoid superellipsoid");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 5.0 0.0 0.0 units box");
    command("set atom 1 mass 1.0");
    command("set atom 2 mass 3.0");
    command("set group all shape 2.0 4.0 6.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("compute ch all chunk/atom type");
    command("compute icc all inertia/chunk ch");
    command("fix h all ave/time 1 1 1 c_icc[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *icc = get_array("icc");
    EXPECT_NEAR(icc[0][0], 10.4, 1.0e-12);
    EXPECT_NEAR(icc[0][1], 26.75, 1.0e-12);
    EXPECT_NEAR(icc[0][2], 22.75, 1.0e-12);
    EXPECT_NEAR(icc[0][3], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][4], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][5], 0.0, 1.0e-12);
}

TEST_F(ComputeInertiaChunkTest, Body)
{
    if (!lammps_config_has_package("BODY")) GTEST_SKIP();

    // single body/nparticle at the origin with a known diagonal inertia
    // tensor (2,3,4); a single chunk must return it unchanged

    const char *datafile = "compute_inertia_chunk_body.data";
    FILE *fp = fopen(datafile, "w");
    ASSERT_NE(fp, nullptr);
    fputs("LAMMPS body nparticle test\n\n"
          "1 atoms\n1 bodies\n1 atom types\n"
          "-10 10 xlo xhi\n-10 10 ylo yhi\n-10 10 zlo zhi\n\n"
          "Atoms\n\n"
          "1 1 1 1.0 0.0 0.0 0.0 0 0 0\n\n"
          "Bodies\n\n"
          "1 1 9\n1\n2.0 3.0 4.0 0.0 0.0 0.0\n0.0 0.0 0.0\n",
          fp);
    fclose(fp);

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style body nparticle 1 1");
    command("read_data " + std::string(datafile));
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("compute ch all chunk/atom type");
    command("compute icc all inertia/chunk ch");
    command("fix h all ave/time 1 1 1 c_icc[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *icc = get_array("icc");
    EXPECT_NEAR(icc[0][0], 2.0, 1.0e-12);
    EXPECT_NEAR(icc[0][1], 3.0, 1.0e-12);
    EXPECT_NEAR(icc[0][2], 4.0, 1.0e-12);
    EXPECT_NEAR(icc[0][3], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][4], 0.0, 1.0e-12);
    EXPECT_NEAR(icc[0][5], 0.0, 1.0e-12);

    remove(datafile);
}

// finite-size particles must contribute their own (spin) angular momentum
// to compute angmom/chunk (issue #3710).

class ComputeAngmomChunkTest : public ComputeChunkTest {
protected:
    void SetUp() override
    {
        testbinary = "ComputeAngmomChunkTest";
        LAMMPSTest::SetUp();
    }
};

TEST_F(ComputeAngmomChunkTest, Sphere)
{
    if (!info->has_style("atom", "sphere")) GTEST_SKIP();

    // single chunk = both spheres (r=1, m=1) spinning at omega=(0,0,5),
    // otherwise at rest -> total angular momentum is pure spin (0,0,4)

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style sphere");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("create_atoms 1 single 2.0 0.0 0.0 units box");
    command("set group all diameter 2.0");
    command("set group all mass 1.0");
    command("set group all omega 0.0 0.0 5.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("compute ch all chunk/atom type");
    command("compute acc all angmom/chunk ch");
    command("fix h all ave/time 1 1 1 c_acc[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *acc = get_array("acc");
    EXPECT_NEAR(acc[0][0], 0.0, 1.0e-12);
    EXPECT_NEAR(acc[0][1], 0.0, 1.0e-12);
    EXPECT_NEAR(acc[0][2], 4.0, 1.0e-12);
}

TEST_F(ComputeAngmomChunkTest, Ellipsoid)
{
    if (!info->has_style("atom", "ellipsoid")) GTEST_SKIP();

    // single ellipsoid (one chunk) at rest with angmom set to (1,2,3)

    BEGIN_HIDE_OUTPUT();
    command("units lj");
    command("atom_style ellipsoid");
    command("boundary f f f");
    command("region box block -20 20 -20 20 -20 20");
    command("create_box 1 box");
    command("create_atoms 1 single 0.0 0.0 0.0 units box");
    command("set group all shape 2.0 4.0 6.0");
    command("set group all mass 1.0");
    command("set group all angmom 1.0 2.0 3.0");
    command("pair_style zero 5.0");
    command("pair_coeff * *");
    command("compute ch all chunk/atom type");
    command("compute acc all angmom/chunk ch");
    command("fix h all ave/time 1 1 1 c_acc[*] mode vector");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    auto *acc = get_array("acc");
    EXPECT_NEAR(acc[0][0], 1.0, 1.0e-12);
    EXPECT_NEAR(acc[0][1], 2.0, 1.0e-12);
    EXPECT_NEAR(acc[0][2], 3.0, 1.0e-12);
}

TEST_F(ComputeChunkTest, ChunkMoleculeWarning)
{
    if (lammps_get_natoms(lmp) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("pair_style lj/cut/coul/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");
    END_HIDE_OUTPUT();

    // compute chunk/atom molecule with a group that excludes *entire*
    // molecules must NOT warn that chunks do not contain all atoms of a
    // molecule, even when an excluded molecule ID is within the range of
    // valid chunk IDs (regression test for issue #5003).  Here molecule 2
    // is excluded entirely, while molecules 1 and 3-6 (max ID = nchunk = 6)
    // are kept intact.  msd/chunk triggers the one-time check at setup.
    // verify both the default (no compress) and the compress code path.
    BEGIN_HIDE_OUTPUT();
    command("group keep molecule 1 3:6");
    command("compute ce1 keep chunk/atom molecule");
    command("compute msd1 all msd/chunk ce1");
    command("compute ce2 keep chunk/atom molecule compress yes");
    command("compute msd2 all msd/chunk ce2");
    END_HIDE_OUTPUT();

    auto output = CAPTURE_OUTPUT([&] { command("run 0 post no"); });
    EXPECT_THAT(output, testing::Not(testing::HasSubstr("do not contain all atoms in molecule")));

    // a molecule that is genuinely *split* between its chunk and the excluded
    // atoms MUST still warn.  excluding atom type 4 keeps most atoms of
    // molecules 1, 2, and 3 but drops their single type-4 atom from the chunk.
    // verify both the default (no compress) and the compress code path.
    BEGIN_HIDE_OUTPUT();
    command("group notype4 type 1 2 3 5");
    command("compute ce3 notype4 chunk/atom molecule");
    command("compute msd3 all msd/chunk ce3");
    command("compute ce4 notype4 chunk/atom molecule compress yes");
    command("compute msd4 all msd/chunk ce4");
    END_HIDE_OUTPUT();

    output = CAPTURE_OUTPUT([&] { command("run 0 post no"); });
    EXPECT_THAT(output, testing::HasSubstr("do not contain all atoms in molecule"));
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
