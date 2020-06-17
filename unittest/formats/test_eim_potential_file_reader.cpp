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

#include "MANYBODY/pair_eim.h"
#include "lammps.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <mpi.h>

using namespace LAMMPS_NS;

class EIMPotentialFileReaderTest : public ::testing::Test {
protected:
    LAMMPS *lmp;
    PairEIM::Setfl setfl;
    static const int nelements = 9;

    void SetUp() override
    {
        const char *args[] = {
            "PotentialFileReaderTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);
        ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        ::testing::internal::GetCapturedStdout();

        int npair        = nelements * (nelements + 1) / 2;
        setfl.ielement   = new int[nelements];
        setfl.mass       = new double[nelements];
        setfl.negativity = new double[nelements];
        setfl.ra         = new double[nelements];
        setfl.ri         = new double[nelements];
        setfl.Ec         = new double[nelements];
        setfl.q0         = new double[nelements];
        setfl.rcutphiA   = new double[npair];
        setfl.rcutphiR   = new double[npair];
        setfl.Eb         = new double[npair];
        setfl.r0         = new double[npair];
        setfl.alpha      = new double[npair];
        setfl.beta       = new double[npair];
        setfl.rcutq      = new double[npair];
        setfl.Asigma     = new double[npair];
        setfl.rq         = new double[npair];
        setfl.rcutsigma  = new double[npair];
        setfl.Ac         = new double[npair];
        setfl.zeta       = new double[npair];
        setfl.rs         = new double[npair];
        setfl.tp         = new int[npair];
    }

    void TearDown() override
    {
        delete[] setfl.ielement;
        delete[] setfl.mass;
        delete[] setfl.negativity;
        delete[] setfl.ra;
        delete[] setfl.ri;
        delete[] setfl.Ec;
        delete[] setfl.q0;
        delete[] setfl.rcutphiA;
        delete[] setfl.rcutphiR;
        delete[] setfl.Eb;
        delete[] setfl.r0;
        delete[] setfl.alpha;
        delete[] setfl.beta;
        delete[] setfl.rcutq;
        delete[] setfl.Asigma;
        delete[] setfl.rq;
        delete[] setfl.rcutsigma;
        delete[] setfl.Ac;
        delete[] setfl.zeta;
        delete[] setfl.rs;
        delete[] setfl.tp;

        ::testing::internal::CaptureStdout();
        delete lmp;
        ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(EIMPotentialFileReaderTest, global_line)
{
    ::testing::internal::CaptureStdout();
    EIMPotentialFileReader reader(lmp, "ffield.eim");
    ::testing::internal::GetCapturedStdout();

    reader.get_global(&setfl);
    ASSERT_DOUBLE_EQ(setfl.division, 2.0);
    ASSERT_DOUBLE_EQ(setfl.rbig, -1.645);
    ASSERT_DOUBLE_EQ(setfl.rsmall, 1.645);
}

TEST_F(EIMPotentialFileReaderTest, element_line_sequential)
{
    ::testing::internal::CaptureStdout();
    EIMPotentialFileReader reader(lmp, "ffield.eim");
    ::testing::internal::GetCapturedStdout();

    reader.get_element(&setfl, 0, "Li");
    ASSERT_EQ(setfl.ielement[0], 3);
    ASSERT_DOUBLE_EQ(setfl.mass[0], 6.9410e+00);
    ASSERT_DOUBLE_EQ(setfl.negativity[0], 9.8000e-01);
    ASSERT_DOUBLE_EQ(setfl.ra[0], 1.1220e+00);
    ASSERT_DOUBLE_EQ(setfl.ri[0], 1.1220e+00);
    ASSERT_DOUBLE_EQ(setfl.Ec[0], -1.6500e+00);
    ASSERT_DOUBLE_EQ(setfl.q0[0], 0.0000e+00);

    reader.get_element(&setfl, 1, "Na");
    ASSERT_EQ(setfl.ielement[1], 11);
    ASSERT_DOUBLE_EQ(setfl.mass[1], 2.2990e+01);
    ASSERT_DOUBLE_EQ(setfl.negativity[1], 9.3000e-01);
    ASSERT_DOUBLE_EQ(setfl.ra[1], 1.3690e+00);
    ASSERT_DOUBLE_EQ(setfl.ri[1], 1.3690e+00);
    ASSERT_DOUBLE_EQ(setfl.Ec[1], -1.1100e+00);
    ASSERT_DOUBLE_EQ(setfl.q0[1], 0.0000e+00);
}

TEST_F(EIMPotentialFileReaderTest, element_line_random)
{
    ::testing::internal::CaptureStdout();
    EIMPotentialFileReader reader(lmp, "ffield.eim");
    ::testing::internal::GetCapturedStdout();

    reader.get_element(&setfl, 0, "Id");
    ASSERT_EQ(setfl.ielement[0], 53);
    ASSERT_DOUBLE_EQ(setfl.mass[0], 1.2690e+02);
    ASSERT_DOUBLE_EQ(setfl.negativity[0], 2.6600e+00);
    ASSERT_DOUBLE_EQ(setfl.ra[0], 1.8500e+00);
    ASSERT_DOUBLE_EQ(setfl.ri[0], 1.8500e+00);
    ASSERT_DOUBLE_EQ(setfl.Ec[0], -1.1100e+00);
    ASSERT_DOUBLE_EQ(setfl.q0[0], 0.0000e+00);
}

TEST_F(EIMPotentialFileReaderTest, pair_line)
{
    ::testing::internal::CaptureStdout();
    EIMPotentialFileReader reader(lmp, "ffield.eim");
    ::testing::internal::GetCapturedStdout();

    reader.get_pair(&setfl, 0, "Li", "Li");
    ASSERT_DOUBLE_EQ(setfl.rcutphiA[0], 6.0490e+00);
    ASSERT_DOUBLE_EQ(setfl.rcutphiR[0], 6.0490e+00);
    ASSERT_DOUBLE_EQ(setfl.Eb[0], -2.5330e-01);
    ASSERT_DOUBLE_EQ(setfl.r0[0], 3.6176e+00);
    ASSERT_DOUBLE_EQ(setfl.alpha[0], 7.5536e+00);
    ASSERT_DOUBLE_EQ(setfl.beta[0], 3.5017e+00);
    ASSERT_DOUBLE_EQ(setfl.rcutq[0], 0.0000e+00);
    ASSERT_DOUBLE_EQ(setfl.Asigma[0], 2.1778e-02);
    ASSERT_DOUBLE_EQ(setfl.rq[0], 2.0000e+00);
    ASSERT_DOUBLE_EQ(setfl.rcutsigma[0], 7.0637e+00);
    ASSERT_DOUBLE_EQ(setfl.Ac[0], 3.3271e-01);
    ASSERT_DOUBLE_EQ(setfl.zeta[0], 6.0000e-01);
    ASSERT_DOUBLE_EQ(setfl.rs[0], 2.0000e+00);
    ASSERT_EQ(setfl.tp[0], 1);
}

TEST_F(EIMPotentialFileReaderTest, pair_identical)
{
    ::testing::internal::CaptureStdout();
    EIMPotentialFileReader reader(lmp, "ffield.eim");
    ::testing::internal::GetCapturedStdout();

    reader.get_pair(&setfl, 0, "Li", "Na");
    reader.get_pair(&setfl, 1, "Na", "Li");
    ASSERT_DOUBLE_EQ(setfl.rcutphiA[0], setfl.rcutphiA[1]);
    ASSERT_DOUBLE_EQ(setfl.rcutphiR[0], setfl.rcutphiR[1]);
    ASSERT_DOUBLE_EQ(setfl.Eb[0], setfl.Eb[1]);
    ASSERT_DOUBLE_EQ(setfl.r0[0], setfl.r0[1]);
    ASSERT_DOUBLE_EQ(setfl.alpha[0], setfl.alpha[1]);
    ASSERT_DOUBLE_EQ(setfl.beta[0], setfl.beta[1]);
    ASSERT_DOUBLE_EQ(setfl.rcutq[0], setfl.rcutq[1]);
    ASSERT_DOUBLE_EQ(setfl.Asigma[0], setfl.Asigma[1]);
    ASSERT_DOUBLE_EQ(setfl.rq[0], setfl.rq[1]);
    ASSERT_DOUBLE_EQ(setfl.rcutsigma[0], setfl.rcutsigma[1]);
    ASSERT_DOUBLE_EQ(setfl.Ac[0], setfl.Ac[1]);
    ASSERT_DOUBLE_EQ(setfl.zeta[0], setfl.zeta[1]);
    ASSERT_DOUBLE_EQ(setfl.rs[0], setfl.rs[1]);
    ASSERT_EQ(setfl.tp[0], setfl.tp[1]);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}
