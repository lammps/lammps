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
#include <cmath>
#include <cstdio>
#include <cstdlib>
/* initialize atoms on fcc lattice in parallel fashion */

#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)


int create_system(System &system, int nx, int ny, int nz, double rho)
{
  /* Box Setup */

  double lattice = pow((4.0 / rho), (1.0 / 3.0));
  system.box.xprd = nx * lattice;
  system.box.yprd = ny * lattice;
  system.box.zprd = nz * lattice;
  system.box.xlo = 0;
  system.box.ylo = 0;
  system.box.zlo = 0;
  system.box.xhi = system.box.xprd;
  system.box.yhi = system.box.yprd;
  system.box.zhi = system.box.zprd;


  int ghost_dist = int(system.neigh_cut/lattice) + 1;

  /* total # of atoms */

  system.nlocal = 4 * nx * ny * nz;
  system.nghost = 4 * (nx + 2 * ghost_dist) *
                      (ny + 2 * ghost_dist) *
                      (nz + 2 * ghost_dist) -
                      system.nlocal;
  system.natoms = system.nlocal + system.nghost;

  system.d_x = t_x_array("X",system.natoms);
  system.h_x = Kokkos::create_mirror_view(system.d_x);
  system.f = t_f_array("F",system.natoms);

  /* determine loop bounds of lattice subsection that overlaps my sub-box
     insure loop bounds do not exceed nx,ny,nz */

  double alat = pow((4.0 / rho), (1.0 / 3.0));
  int ilo = static_cast<int>(system.box.xlo / (0.5 * alat) - 1);
  int ihi = static_cast<int>(system.box.xhi / (0.5 * alat) + 1);
  int jlo = static_cast<int>(system.box.ylo / (0.5 * alat) - 1);
  int jhi = static_cast<int>(system.box.yhi / (0.5 * alat) + 1);
  int klo = static_cast<int>(system.box.zlo / (0.5 * alat) - 1);
  int khi = static_cast<int>(system.box.zhi / (0.5 * alat) + 1);

  ilo = MAX(ilo, 0);
  ihi = MIN(ihi, 2 * nx - 1);
  jlo = MAX(jlo, 0);
  jhi = MIN(jhi, 2 * ny - 1);
  klo = MAX(klo, 0);
  khi = MIN(khi, 2 * nz - 1);



  /* generates positions of atoms on fcc sublattice*/

  srand(3718273);
  /* create non-ghost atoms */
  {
    double xtmp, ytmp, ztmp;
    int sx = 0;
    int sy = 0;
    int sz = 0;
    int ox = 0;
    int oy = 0;
    int oz = 0;
    int subboxdim = 8;

    int n = 0;
    int iflag = 0;

    while(oz * subboxdim <= khi) {
      const int k = oz * subboxdim + sz;
      const int j = oy * subboxdim + sy;
      const int i = ox * subboxdim + sx;

      if(iflag) continue;

      if(((i + j + k) % 2 == 0) &&
          (i >= ilo) && (i <= ihi) &&
          (j >= jlo) && (j <= jhi) &&
          (k >= klo) && (k <= khi)) {

        const int nold = n;
        while(nold == n) {
          xtmp = 0.5 * alat * i + system.delta/1000*(rand()%1000-500);
          ytmp = 0.5 * alat * j + system.delta/1000*(rand()%1000-500);
          ztmp = 0.5 * alat * k + system.delta/1000*(rand()%1000-500);

          if(xtmp >= system.box.xlo && xtmp < system.box.xhi &&
              ytmp >= system.box.ylo && ytmp < system.box.yhi &&
              ztmp >= system.box.zlo && ztmp < system.box.zhi) {
            system.h_x(n,0) = xtmp;
            system.h_x(n,1) = ytmp;
            system.h_x(n,2) = ztmp;
            n++;
          }
        }
      }

      sx++;

      if(sx == subboxdim) {
        sx = 0;
        sy++;
      }

      if(sy == subboxdim) {
        sy = 0;
        sz++;
      }

      if(sz == subboxdim) {
        sz = 0;
        ox++;
      }

      if(ox * subboxdim > ihi) {
        ox = 0;
        oy++;
      }

      if(oy * subboxdim > jhi) {
        oy = 0;
        oz++;
      }
    }

    /* check that correct # of atoms were created */

    if(system.nlocal != n) {
      printf("Created incorrect # of atoms\n");

      return 1;
    }
  }

  /* create ghost atoms */

  {
    double xtmp, ytmp, ztmp;

    int ilo_g = ilo - 2 * ghost_dist;
    int jlo_g = jlo - 2 * ghost_dist;
    int klo_g = klo - 2 * ghost_dist;
    int ihi_g = ihi + 2 * ghost_dist;
    int jhi_g = jhi + 2 * ghost_dist;
    int khi_g = khi + 2 * ghost_dist;

    int subboxdim = 8;
    int sx = 0;
    int sy = 0;
    int sz = 0;
    int ox = subboxdim * ilo_g;
    int oy = subboxdim * jlo_g;
    int oz = subboxdim * klo_g;

    int n = system.nlocal;
    int iflag = 0;


    while(oz * subboxdim <= khi_g) {
      const int k = oz * subboxdim + sz;
      const int j = oy * subboxdim + sy;
      const int i = ox * subboxdim + sx;

      if(iflag) continue;

      if(((i + j + k) % 2 == 0) &&
          (i >= ilo_g) && (i <= ihi_g) &&
          (j >= jlo_g) && (j <= jhi_g) &&
          (k >= klo_g) && (k <= khi_g) &&
          ((i < ilo) || (i > ihi) ||
           (j < jlo) || (j > jhi) ||
           (k < klo) || (k > khi))
          ) {

        xtmp = 0.5 * alat * i;
        ytmp = 0.5 * alat * j;
        ztmp = 0.5 * alat * k;

        system.h_x(n,0) = xtmp + system.delta/1000*(rand()%1000-500);;
        system.h_x(n,1) = ytmp + system.delta/1000*(rand()%1000-500);;
        system.h_x(n,2) = ztmp + system.delta/1000*(rand()%1000-500);;
        n++;
      }

      sx++;

      if(sx == subboxdim) {
        sx = 0;
        sy++;
      }

      if(sy == subboxdim) {
        sy = 0;
        sz++;
      }

      if(sz == subboxdim) {
        sz = 0;
        ox++;
        //printf("%i %i %i // %i %i %i\n",ox,oy,oz,i,j,k);
      }

      if(ox * subboxdim > ihi_g) {
        ox = subboxdim * ilo_g;
        oy++;
      }

      if(oy * subboxdim > jhi_g) {
        oy = subboxdim * jlo_g;
        oz++;
      }
    }
  }

  Kokkos::deep_copy(system.d_x,system.h_x);
  return 0;
}

