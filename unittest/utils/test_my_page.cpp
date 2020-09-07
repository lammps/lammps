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

#include "lmptype.h"
#include "my_page.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace LAMMPS_NS;

TEST(MyPage, int_default) {
    MyPage<int> p;

    // default init. maxchunk=1, pagesize=1024
    int rv = p.init();
    ASSERT_EQ(rv,0);

    ASSERT_EQ(p.ndatum,0);
    ASSERT_EQ(p.nchunk,0);

    int *iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr,p.vget());
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(0,p.status());
    ASSERT_EQ(p.ndatum,1);
    ASSERT_EQ(p.nchunk,1);
    ASSERT_EQ(iptr,p.vget());
    // use too large chunk size
    p.vgot(2);
    ASSERT_EQ(1,p.status());

    p.reset();
    ASSERT_EQ(p.ndatum,0);
    ASSERT_EQ(p.nchunk,0);

    iptr = p.vget();
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(iptr,p.get());
    ++iptr;
    ASSERT_EQ(iptr,p.get(1));
    ASSERT_EQ(p.ndatum,3);
    ASSERT_EQ(p.nchunk,3);
}

TEST(MyPage, int_custom) {
    MyPage<int> p;

    // default init. maxchunk=16, pagesize=256
    int rv = p.init(16,256);
    ASSERT_EQ(rv,0);

    ASSERT_EQ(p.ndatum,0);
    ASSERT_EQ(p.nchunk,0);

    int *iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr,p.vget());
    p.vgot(16);
    iptr += 16;
    ASSERT_EQ(0,p.status());
    ASSERT_EQ(p.ndatum,16);
    ASSERT_EQ(p.nchunk,1);

    // use too large chunk size
    ASSERT_EQ(iptr,p.vget());
    p.vgot(32);
    ASSERT_EQ(1,p.status());

    p.reset();
    ASSERT_EQ(p.ndatum,0);
    ASSERT_EQ(p.nchunk,0);

    iptr = p.vget();
    p.vgot(16);
    iptr = p.vget();
    p.vgot(4);
    iptr += 4;
    ASSERT_EQ(iptr,p.get());
    ++iptr;
    ASSERT_EQ(iptr,p.get(16));
    ASSERT_EQ(p.ndatum,37);
    ASSERT_EQ(p.nchunk,4);
}
