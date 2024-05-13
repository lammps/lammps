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

#include "lmptype.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

using namespace LAMMPS_NS;

TEST(Types, ubuf)
{
    double buf[3];
    double d1 = 0.1;
    int i1 = -10;
#if defined(LAMMPS_SMALLSMALL)
    bigint b1 = 2048;
#else
    bigint b1 = (1L << 58) + (1L << 50);
#endif
    buf[0] = d1;
    buf[1] = ubuf(i1).d;
    buf[2] = ubuf(b1).d;

    EXPECT_EQ(d1, buf[0]);
    EXPECT_EQ(i1, (int)ubuf(buf[1]).i);
    EXPECT_EQ(b1, (bigint)ubuf(buf[2]).i);
}

TEST(Types, multitype)
{
    multitype m[7];
    int64_t b1 = (3L << 48) - 1;
    int i1    = 20;
    double d1 = 0.1;

    m[0] = b1;
    m[1] = i1;
    m[2] = d1;

    m[3] = (bigint) -((1L << 40) + (1L << 50));
    m[4] = -1023;
    m[5] = -2.225;

    EXPECT_EQ(m[0].type, multitype::LAMMPS_INT64);
    EXPECT_EQ(m[1].type, multitype::LAMMPS_INT);
    EXPECT_EQ(m[2].type, multitype::LAMMPS_DOUBLE);

#if defined(LAMMPS_SMALLSMALL)
    EXPECT_EQ(m[3].type, multitype::LAMMPS_INT);
#else
    EXPECT_EQ(m[3].type, multitype::LAMMPS_INT64);
#endif
    EXPECT_EQ(m[4].type, multitype::LAMMPS_INT);
    EXPECT_EQ(m[5].type, multitype::LAMMPS_DOUBLE);
    EXPECT_EQ(m[6].type, multitype::LAMMPS_NONE);

    EXPECT_EQ(m[0].data.b, b1);
    EXPECT_EQ(m[1].data.i, i1);
    EXPECT_EQ(m[2].data.d, d1);

#if !defined(LAMMPS_SMALLSMALL)
    EXPECT_EQ(m[3].data.b,  -((1L << 40) + (1L << 50)));
#endif
    EXPECT_EQ(m[4].data.i, -1023);
    EXPECT_EQ(m[5].data.d, -2.225);
}
