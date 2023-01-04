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
#include "my_page.h"
#include "my_pool_chunk.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace LAMMPS_NS;

TEST(MyPage, int)
{
    MyPage<int> p;

    // default init. maxchunk=1, pagesize=1024
    int rv = p.init();
    ASSERT_EQ(rv, 0);

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    int *iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr, p.vget());
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 1);
    ASSERT_EQ(p.nchunk, 1);
    ASSERT_EQ(iptr, p.vget());
    // use too large chunk size
    p.vgot(2);
    ASSERT_EQ(1, p.status());

    p.reset();
    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(iptr, p.get());
    ++iptr;
    ASSERT_EQ(iptr, p.get(1));
    ASSERT_EQ(p.ndatum, 3);
    ASSERT_EQ(p.nchunk, 3);

    // restart with custom init. maxchunk=16, pagesize=256
    rv = p.init(16, 64, 2);
    ASSERT_EQ(rv, 0);

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr, p.vget());
    p.vgot(16);
    iptr += 16;
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 16);
    ASSERT_EQ(p.nchunk, 1);

    // use too large chunk size
    ASSERT_EQ(iptr, p.vget());
    p.vgot(32);
    ASSERT_EQ(1, p.status());

    p.reset();
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    p.vgot(16);
    iptr = p.vget();
    p.vgot(4);
    iptr += 4;
    ASSERT_EQ(iptr, p.get());
    ++iptr;
    ASSERT_EQ(iptr, p.get(16));
    ASSERT_DOUBLE_EQ(p.size(), (double)sizeof(int) * 128.0);
    ASSERT_EQ(p.ndatum, 37);
    ASSERT_EQ(p.nchunk, 4);
    p.get(16);
    p.get(16);
    // allocation on the same page
    iptr = p.get(16);
    iptr += 16;
    ASSERT_EQ(iptr, p.get(16));
    // allocation on different pages
    p.get(16);
    iptr += 16;
    ASSERT_NE(iptr, p.get(16));
    ASSERT_DOUBLE_EQ(p.size(), (double)sizeof(int) * 256.0);
    ASSERT_EQ(p.ndatum, 133);
    ASSERT_EQ(p.nchunk, 10);
}

TEST(MyPage, double)
{
    MyPage<double> p;

    // default init. maxchunk=1, pagesize=1024
    int rv = p.init();
    ASSERT_EQ(rv, 0);

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    double *iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr, p.vget());
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 1);
    ASSERT_EQ(p.nchunk, 1);
    ASSERT_EQ(iptr, p.vget());
    // use too large chunk size
    p.vgot(2);
    ASSERT_EQ(1, p.status());

    p.reset();
    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(iptr, p.get());
    ++iptr;
    ASSERT_EQ(iptr, p.get(1));
    ASSERT_EQ(p.ndatum, 3);
    ASSERT_EQ(p.nchunk, 3);

    // restart with custom init. maxchunk=16, pagesize=256
    rv = p.init(16, 64, 2);
    ASSERT_EQ(rv, 0);

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr, p.vget());
    p.vgot(16);
    iptr += 16;
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 16);
    ASSERT_EQ(p.nchunk, 1);

    // use too large chunk size
    ASSERT_EQ(iptr, p.vget());
    p.vgot(32);
    ASSERT_EQ(1, p.status());

    p.reset();
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    p.vgot(16);
    iptr = p.vget();
    p.vgot(4);
    iptr += 4;
    ASSERT_EQ(iptr, p.get());
    ++iptr;
    ASSERT_EQ(iptr, p.get(16));
    ASSERT_DOUBLE_EQ(p.size(), (double)sizeof(double) * 128.0);
    ASSERT_EQ(p.ndatum, 37);
    ASSERT_EQ(p.nchunk, 4);
    p.get(16);
    p.get(16);
    // allocation on the same page
    iptr = p.get(16);
    iptr += 16;
    ASSERT_EQ(iptr, p.get(16));
    // allocation on different pages
    p.get(16);
    iptr += 16;
    ASSERT_NE(iptr, p.get(16));
    ASSERT_DOUBLE_EQ(p.size(), (double)sizeof(double) * 256.0);
    ASSERT_EQ(p.ndatum, 133);
    ASSERT_EQ(p.nchunk, 10);
}

TEST(MyPage, bigint)
{
    MyPage<bigint> p;

    // default init. maxchunk=1, pagesize=1024
    int rv = p.init();
    ASSERT_EQ(rv, 0);

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    bigint *iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr, p.vget());
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 1);
    ASSERT_EQ(p.nchunk, 1);
    ASSERT_EQ(iptr, p.vget());
    // use too large chunk size
    p.vgot(2);
    ASSERT_EQ(1, p.status());

    p.reset();
    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    p.vgot(1);
    ++iptr;
    ASSERT_EQ(iptr, p.get());
    ++iptr;
    ASSERT_EQ(iptr, p.get(1));
    ASSERT_EQ(p.ndatum, 3);
    ASSERT_EQ(p.nchunk, 3);

    // restart with custom init. maxchunk=16, pagesize=256
    rv = p.init(16, 64, 2);
    ASSERT_EQ(rv, 0);

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    // second call to vget() should give same pointer without vgot()
    ASSERT_EQ(iptr, p.vget());
    p.vgot(16);
    iptr += 16;
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 16);
    ASSERT_EQ(p.nchunk, 1);

    // use too large chunk size
    ASSERT_EQ(iptr, p.vget());
    p.vgot(32);
    ASSERT_EQ(1, p.status());

    p.reset();
    ASSERT_EQ(0, p.status());
    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);

    iptr = p.vget();
    p.vgot(16);
    iptr = p.vget();
    p.vgot(4);
    iptr += 4;
    ASSERT_EQ(iptr, p.get());
    ++iptr;
    ASSERT_EQ(iptr, p.get(16));
    ASSERT_DOUBLE_EQ(p.size(), (double)sizeof(bigint) * 128.0);
    ASSERT_EQ(p.ndatum, 37);
    ASSERT_EQ(p.nchunk, 4);
    p.get(16);
    p.get(16);
    // allocation on the same page
    iptr = p.get(16);
    iptr += 16;
    ASSERT_EQ(iptr, p.get(16));
    // allocation on different pages
    p.get(16);
    iptr += 16;
    ASSERT_NE(iptr, p.get(16));
    ASSERT_DOUBLE_EQ(p.size(), (double)sizeof(bigint) * 256.0);
    ASSERT_EQ(p.ndatum, 133);
    ASSERT_EQ(p.nchunk, 10);
}

TEST(MyPoolChunk, int)
{
    // defaults to minchunk=1, maxchunk=1, nbin=1, chunksperpage=1024, pagedelta=1
    MyPoolChunk<int> p;

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);
    ASSERT_EQ(p.size(), 0.0);

    int idx   = ~0x0000;
    int *iptr = p.get(idx);
    ASSERT_NE(iptr, nullptr);
    ASSERT_EQ(idx, 0);

    iptr = p.get(1, idx);
    ASSERT_NE(iptr, nullptr);
    ASSERT_EQ(idx, 1);
    // we have only one page allocated
    ASSERT_EQ(p.size(), 1024 * sizeof(int) + 1024 * sizeof(int) + sizeof(void *) + sizeof(int));
    ASSERT_EQ(p.ndatum, 2);
    ASSERT_EQ(p.nchunk, 2);

    p.put(0);
    ASSERT_EQ(p.ndatum, 1);
    ASSERT_EQ(p.nchunk, 1);

    iptr = p.get(2, idx);
    ASSERT_EQ(iptr, nullptr);
    ASSERT_EQ(p.status(), 3);
    ASSERT_EQ(p.ndatum, 1);
    ASSERT_EQ(p.nchunk, 1);
}

TEST(MyPoolChunk, double)
{
    // defaults to minchunk=1, maxchunk=1, nbin=1, chunksperpage=1024, pagedelta=1
    MyPoolChunk<double> p;

    ASSERT_EQ(p.ndatum, 0);
    ASSERT_EQ(p.nchunk, 0);
    ASSERT_EQ(p.size(), 0.0);

    int idx      = ~0x0000;
    double *dptr = p.get(idx);
    ASSERT_NE(dptr, nullptr);
    ASSERT_EQ(idx, 0);

    dptr = p.get(1, idx);
    ASSERT_NE(dptr, nullptr);
    ASSERT_EQ(idx, 1);
    // we have only one page allocated
    ASSERT_EQ(p.size(), 1024 * sizeof(int) + 1024 * sizeof(double) + sizeof(void *) + sizeof(int));

    p.put(0);
    ASSERT_EQ(p.ndatum, 1);
    ASSERT_EQ(p.nchunk, 1);

    dptr = p.get(2, idx);
    ASSERT_EQ(dptr, nullptr);
    ASSERT_EQ(p.status(), 3);
    ASSERT_EQ(p.ndatum, 1);
    ASSERT_EQ(p.nchunk, 1);
}
