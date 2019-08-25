/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include <system.h>
#include <cstdio>
#include <Kokkos_Core.hpp>

#define SMALL 1.0e-6
#define FACTOR 0.999

/* BinningFunctor puts atoms into bins of the simulation box
 * Neighborlists are then created by checking only distances of atoms
 * in adjacent bins. That makes neighborlist construction a O(N) operation.
 */

struct BinningFunctor {
  typedef t_int_2d::execution_space execution_space;

  System s;

  int atoms_per_bin;

  BinningFunctor(System _s): s(_s) {
    atoms_per_bin = s.bins.dimension_1();
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const
  {
    const int ibin = coord2bin(s.d_x(i, 0), s.d_x(i, 1), s.d_x(i, 2));

    const int ac = Kokkos::atomic_fetch_add(&s.bincount[ibin], 1);

    if(ac < atoms_per_bin) {
      s.bins(ibin, ac) = i;
    } else if(s.d_resize(0) < ac) {
      s.d_resize(0) = ac;
    }
  }

  KOKKOS_INLINE_FUNCTION
  int coord2bin(double x, double y, double z) const
  {
    int ix, iy, iz;

    if(x >= s.box.xprd)
      ix = (int)((x - s.box.xprd) * s.bininvx) + s.nbinx - s.mbinxlo;
    else if(x >= 0.0)
      ix = (int)(x * s.bininvx) - s.mbinxlo;
    else
      ix = (int)(x * s.bininvx) - s.mbinxlo - 1;

    if(y >= s.box.yprd)
      iy = (int)((y - s.box.yprd) * s.bininvy) + s.nbiny - s.mbinylo;
    else if(y >= 0.0)
      iy = (int)(y * s.bininvy) - s.mbinylo;
    else
      iy = (int)(y * s.bininvy) - s.mbinylo - 1;

    if(z >= s.box.zprd)
      iz = (int)((z - s.box.zprd) * s.bininvz) + s.nbinz - s.mbinzlo;
    else if(z >= 0.0)
      iz = (int)(z * s.bininvz) - s.mbinzlo;
    else
      iz = (int)(z * s.bininvz) - s.mbinzlo - 1;

    return (iz * s.mbiny * s.mbinx + iy * s.mbinx + ix + 1);
  }
};

/* Build the actual neighborlist*/

struct BuildFunctor {

  typedef t_int_2d::execution_space execution_space;

  System s;

  int maxneighs;
  BuildFunctor(System _s): s(_s) {
    maxneighs = s.neighbors.dimension_1();
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const
  {
    int n = 0;

    const t_int_1d_const_um bincount_c = s.bincount;

    const double xtmp = s.d_x(i, 0);
    const double ytmp = s.d_x(i, 1);
    const double ztmp = s.d_x(i, 2);

    const int ibin = coord2bin(xtmp, ytmp, ztmp);

    // loop over all bins in neighborhood (includes ibin)
    for(int k = 0; k < s.nstencil; k++) {
      const int jbin = ibin + s.d_stencil[k];

      // get subview of jbin
      const t_int_1d_const_um loc_bin =
          Kokkos::subview(s.bins,jbin,Kokkos::ALL());

      if(ibin == jbin)
        for(int m = 0; m < bincount_c[jbin]; m++) {
          const int j = loc_bin[m];

          //for same bin as atom i skip j if i==j
          if (j == i) continue;

          const double delx = xtmp - s.d_x(j, 0);
          const double dely = ytmp - s.d_x(j, 1);
          const double delz = ztmp - s.d_x(j, 2);
          const double rsq = delx * delx + dely * dely + delz * delz;

          if(rsq <= s.neigh_cutsq && n<maxneighs) s.neighbors(i,n++) = j;
        }
      else {
        for(int m = 0; m < bincount_c[jbin]; m++) {
          const int j = loc_bin[m];

          const double delx = xtmp - s.d_x(j, 0);
          const double dely = ytmp - s.d_x(j, 1);
          const double delz = ztmp - s.d_x(j, 2);
          const double rsq = delx * delx + dely * dely + delz * delz;

          if(rsq <= s.neigh_cutsq && n<maxneighs) s.neighbors(i,n++) = j;
        }
      }
    }

    s.numneigh[i] = n;

    if(n >= maxneighs) {
      if(n >= s.d_resize(0)) s.d_resize(0) = n;
    }
  }

  KOKKOS_INLINE_FUNCTION
  int coord2bin(double x, double y, double z) const
  {
    int ix, iy, iz;

    if(x >= s.box.xprd)
      ix = (int)((x - s.box.xprd) * s.bininvx) + s.nbinx - s.mbinxlo;
    else if(x >= 0.0)
      ix = (int)(x * s.bininvx) - s.mbinxlo;
    else
      ix = (int)(x * s.bininvx) - s.mbinxlo - 1;

    if(y >= s.box.yprd)
      iy = (int)((y - s.box.yprd) * s.bininvy) + s.nbiny - s.mbinylo;
    else if(y >= 0.0)
      iy = (int)(y * s.bininvy) - s.mbinylo;
    else
      iy = (int)(y * s.bininvy) - s.mbinylo - 1;

    if(z >= s.box.zprd)
      iz = (int)((z - s.box.zprd) * s.bininvz) + s.nbinz - s.mbinzlo;
    else if(z >= 0.0)
      iz = (int)(z * s.bininvz) - s.mbinzlo;
    else
      iz = (int)(z * s.bininvz) - s.mbinzlo - 1;

    return (iz * s.mbiny * s.mbinx + iy * s.mbinx + ix + 1);
  }
};

/* Reset an array to zero */

struct MemsetZeroFunctor {
  typedef t_x_array::execution_space  execution_space ;
  void* ptr;
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    ((int*)ptr)[i] = 0;
  }
};

/* Calculate distance of two bins */

double bindist(System &s, int i, int j, int k)
{
  double delx, dely, delz;

  if(i > 0)
    delx = (i - 1) * s.binsizex;
  else if(i == 0)
    delx = 0.0;
  else
    delx = (i + 1) * s.binsizex;

  if(j > 0)
    dely = (j - 1) * s.binsizey;
  else if(j == 0)
    dely = 0.0;
  else
    dely = (j + 1) * s.binsizey;

  if(k > 0)
    delz = (k - 1) * s.binsizez;
  else if(k == 0)
    delz = 0.0;
  else
    delz = (k + 1) * s.binsizez;

  return (delx * delx + dely * dely + delz * delz);
}

/* Setup the neighborlist construction
 * Determine binsizes, a stencil for defining adjacency, etc.
 */

void neigh_setup(System &s) {

  s.neigh_cutsq = s.neigh_cut * s.neigh_cut;

  /*
  c bins must evenly divide into box size,
  c   becoming larger than cutneigh if necessary
  c binsize = 1/2 of cutoff is near optimal

  if (flag == 0) {
    nbinx = 2.0 * xprd / cutneigh;
    nbiny = 2.0 * yprd / cutneigh;
    nbinz = 2.0 * zprd / cutneigh;
    if (nbinx == 0) nbinx = 1;
    if (nbiny == 0) nbiny = 1;
    if (nbinz == 0) nbinz = 1;
  }
  */

  s.binsizex = s.box.xprd / s.nbinx;
  s.binsizey = s.box.yprd / s.nbiny;
  s.binsizez = s.box.zprd / s.nbinz;
  s.bininvx = 1.0 / s.binsizex;
  s.bininvy = 1.0 / s.binsizey;
  s.bininvz = 1.0 / s.binsizez;

  double coord = s.box.xlo - s.neigh_cut - SMALL * s.box.xprd;
  s.mbinxlo = static_cast<int>(coord * s.bininvx);

  if(coord < 0.0) s.mbinxlo = s.mbinxlo - 1;

  coord = s.box.xhi + s.neigh_cut + SMALL * s.box.xprd;
  int mbinxhi = static_cast<int>(coord * s.bininvx);

  coord = s.box.ylo - s.neigh_cut - SMALL * s.box.yprd;
  s.mbinylo = static_cast<int>(coord * s.bininvy);

  if(coord < 0.0) s.mbinylo = s.mbinylo - 1;

  coord = s.box.yhi + s.neigh_cut + SMALL * s.box.yprd;
  int mbinyhi = static_cast<int>(coord * s.bininvy);

  coord = s.box.zlo - s.neigh_cut - SMALL * s.box.zprd;
  s.mbinzlo = static_cast<int>(coord * s.bininvz);

  if(coord < 0.0) s.mbinzlo = s.mbinzlo - 1;

  coord = s.box.zhi + s.neigh_cut + SMALL * s.box.zprd;
  int mbinzhi = static_cast<int>(coord * s.bininvz);

  /* extend bins by 1 in each direction to insure stencil coverage */

  s.mbinxlo = s.mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  s.mbinx = mbinxhi - s.mbinxlo + 1;

  s.mbinylo = s.mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  s.mbiny = mbinyhi - s.mbinylo + 1;

  s.mbinzlo = s.mbinzlo - 1;
  mbinzhi = mbinzhi + 1;
  s.mbinz = mbinzhi - s.mbinzlo + 1;

  /*
  compute bin stencil of all bins whose closest corner to central bin
  is within neighbor cutoff
  for partial Newton (newton = 0),
  stencil is all surrounding bins including self
  for full Newton (newton = 1),
  stencil is bins to the "upper right" of central bin, does NOT include self
  next(xyz) = how far the stencil could possibly extend
  factor < 1.0 for special case of LJ benchmark so code will create
  correct-size stencil when there are 3 bins for every 5 lattice spacings
  */

  int nextx = static_cast<int>(s.neigh_cut * s.bininvx);

  if(nextx * s.binsizex < FACTOR * s.neigh_cut) nextx++;

  int nexty = static_cast<int>(s.neigh_cut * s.bininvy);

  if(nexty * s.binsizey < FACTOR * s.neigh_cut) nexty++;

  int nextz = static_cast<int>(s.neigh_cut * s.bininvz);

  if(nextz * s.binsizez < FACTOR * s.neigh_cut) nextz++;

  int nmax = (2 * nextz + 1) * (2 * nexty + 1) * (2 * nextx + 1);
  s.d_stencil = t_int_1d("stencil", nmax);
  s.h_stencil = Kokkos::create_mirror_view(s.d_stencil);
  s.nstencil = 0;
  int kstart = -nextz;

  for(int k = kstart; k <= nextz; k++) {
    for(int j = -nexty; j <= nexty; j++) {
      for(int i = -nextx; i <= nextx; i++) {
        if(bindist(s,i, j, k) < s.neigh_cutsq) {
          s.h_stencil(s.nstencil++) = k * s.mbiny * s.mbinx + j * s.mbinx + i;
        }
      }
    }
  }

  /* Allocate neighbor arrays */

  Kokkos::deep_copy(s.d_stencil, s.h_stencil);
  s.mbins = s.mbinx * s.mbiny * s.mbinz;
  s.bincount = t_int_1d("bincount", s.mbins);
  s.bins = t_int_2d("bins", s.mbins, 8);

  s.neighbors = t_neighbors("neighbors",s.natoms,80);
  s.numneigh = t_int_1d("numneigh",s.natoms);
  s.d_resize = t_int_scalar("resize");
  s.h_resize = Kokkos::create_mirror_view(s.d_resize);
}


/* Build the neighborlist
 * This is a try and rerun algorithm for handling the case where the bins array
 * and the neighbors array are not big enough. So if one is too small, it will
 * reallocate and rerun the binnind algorithm or the neighborlist construction.
 */

void neigh_build(System &s) {

  /* Binning of atoms */

  s.h_resize(0) = 1;

  while(s.h_resize(0) > 0) {
    s.h_resize(0) = 0;
    Kokkos::deep_copy(s.d_resize, s.h_resize);

    MemsetZeroFunctor f_zero;
    f_zero.ptr = (void*) s.bincount.ptr_on_device();
    Kokkos::parallel_for(s.mbins, f_zero);
    execution_space().fence();

    BinningFunctor f(s);
    Kokkos::parallel_for(s.natoms, f);
    execution_space().fence();

    /* Check if bins was large enough, if nor reallocated and rerun */

    deep_copy(s.h_resize, s.d_resize);

    if(s.h_resize(0)) {
      int atoms_per_bin = s.h_resize(0)+2;
      s.bins = t_int_2d("bins", s.mbins, atoms_per_bin);
    }
  }

  /* Neighborlist construction */

  s.h_resize(0) = 1;

  while(s.h_resize(0)) {
    s.h_resize(0) = 0;

    Kokkos::deep_copy(s.d_resize, s.h_resize);

    BuildFunctor f(s);
    Kokkos::parallel_for(s.nlocal, f);

    execution_space().fence();

    /* Check if neighbors was large enough, if nor reallocated and rerun */

    deep_copy(s.h_resize, s.d_resize);

    if(s.h_resize(0)) {
      int maxneighs = s.h_resize(0) * 1.2;
      s.neighbors = t_neighbors("neighbors", s.natoms, maxneighs);
    }
  }
}
