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

#include "fix_minimize_kokkos.h"
#include "atom_kokkos.h"
#include "domain.h"
#include "memory_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMinimizeKokkos::FixMinimizeKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixMinimize(lmp, narg, arg)
{
  atomKK = (AtomKokkos *) atom;
}

/* ---------------------------------------------------------------------- */

FixMinimizeKokkos::~FixMinimizeKokkos()
{
  memoryKK->destroy_kokkos(k_vectors,vectors);
  vectors = NULL;
}

/* ----------------------------------------------------------------------
   allocate/initialize memory for a new vector with 3 elements per atom
------------------------------------------------------------------------- */

void FixMinimizeKokkos::add_vector_kokkos()
{
  int n = 3;

  memory->grow(peratom,nvector+1,"minimize:peratom");
  peratom[nvector] = n;

  // d_vectors needs to be LayoutRight for subviews

  k_vectors.sync<LMPDeviceType>();

  memoryKK->grow_kokkos(k_vectors,vectors,nvector+1,atom->nmax*n,
                      "minimize:vectors");
  d_vectors = k_vectors.d_view;
  h_vectors = k_vectors.h_view;

  k_vectors.modify<LMPDeviceType>();

  nvector++;
}

/* ----------------------------------------------------------------------
   return a pointer to the Mth vector
------------------------------------------------------------------------- */

DAT::t_ffloat_1d FixMinimizeKokkos::request_vector_kokkos(int m)
{
  k_vectors.sync<LMPDeviceType>();

  return Kokkos::subview(d_vectors,m,Kokkos::ALL);
}

/* ----------------------------------------------------------------------
   reset x0 for atoms that moved across PBC via reneighboring in line search
   x0 = 1st vector
   must do minimum_image using original box stored at beginning of line search
   swap & set_global_box() change to original box, then restore current box
------------------------------------------------------------------------- */

void FixMinimizeKokkos::reset_coords()
{
  box_swap();
  domain->set_global_box();

  int nlocal = atom->nlocal;

  atomKK->sync(Device,X_MASK);
  k_vectors.sync<LMPDeviceType>();

  {
    // local variables for lambda capture

    auto triclinic = domain->triclinic;
    auto xperiodic = domain->xperiodic;
    auto xprd_half = domain->xprd_half;
    auto xprd = domain->xprd;
    auto yperiodic = domain->yperiodic;
    auto yprd_half = domain->yprd_half;
    auto yprd = domain->yprd;
    auto zperiodic = domain->zperiodic;
    auto zprd_half = domain->zprd_half;
    auto zprd = domain->zprd;
    auto xy = domain->xy;
    auto xz = domain->xz;
    auto yz = domain->yz;
    auto l_x = atomKK->k_x.d_view;
    auto l_x0 = Kokkos::subview(d_vectors,0,Kokkos::ALL);

    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(const int& i) {
      const int n = i*3;
      double dx0 = l_x(i,0) - l_x0[n];
      double dy0 = l_x(i,1) - l_x0[n+1];
      double dz0 = l_x(i,2) - l_x0[n+2];
      double dx = dx0;
      double dy = dy0;
      double dz = dz0;
      // domain->minimum_image(dx,dy,dz);
      {
        if (triclinic == 0) {
          if (xperiodic) {
            if (fabs(dx) > xprd_half) {
              if (dx < 0.0) dx += xprd;
              else dx -= xprd;
            }
          }
          if (yperiodic) {
            if (fabs(dy) > yprd_half) {
              if (dy < 0.0) dy += yprd;
              else dy -= yprd;
            }
          }
          if (zperiodic) {
            if (fabs(dz) > zprd_half) {
              if (dz < 0.0) dz += zprd;
              else dz -= zprd;
            }
          }

        } else {
          if (zperiodic) {
            if (fabs(dz) > zprd_half) {
              if (dz < 0.0) {
                dz += zprd;
                dy += yz;
                dx += xz;
              } else {
                dz -= zprd;
                dy -= yz;
                dx -= xz;
              }
            }
          }
          if (yperiodic) {
            if (fabs(dy) > yprd_half) {
              if (dy < 0.0) {
                dy += yprd;
                dx += xy;
              } else {
                dy -= yprd;
                dx -= xy;
              }
            }
          }
          if (xperiodic) {
            if (fabs(dx) > xprd_half) {
              if (dx < 0.0) dx += xprd;
              else dx -= xprd;
            }
          }
        }
      } // end domain->minimum_image(dx,dy,dz);
      if (dx != dx0) l_x0[n] = l_x(i,0) - dx;
      if (dy != dy0) l_x0[n+1] = l_x(i,1) - dy;
      if (dz != dz0) l_x0[n+2] = l_x(i,2) - dz;
    });
  }
  k_vectors.modify<LMPDeviceType>();

  box_swap();
  domain->set_global_box();
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixMinimizeKokkos::grow_arrays(int nmax)
{
  k_vectors.sync<LMPDeviceType>();
  memoryKK->grow_kokkos(k_vectors,vectors,nvector,3*nmax,"minimize:vector");
  d_vectors = k_vectors.d_view;
  h_vectors = k_vectors.h_view;
  k_vectors.modify<LMPDeviceType>();
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixMinimizeKokkos::copy_arrays(int i, int j, int /*delflag*/)
{
  int m,iper,nper,ni,nj;

  k_vectors.sync<LMPHostType>();

  for (m = 0; m < nvector; m++) {
    nper = 3;
    ni = nper*i;
    nj = nper*j;
    for (iper = 0; iper < nper; iper++) h_vectors(m,nj++) = h_vectors(m,ni++);
  }

  k_vectors.modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixMinimizeKokkos::pack_exchange(int i, double *buf)
{
  int m,iper,nper,ni;

  k_vectors.sync<LMPHostType>();

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*i;
    for (iper = 0; iper < nper; iper++) buf[n++] = h_vectors(m,ni++);
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixMinimizeKokkos::unpack_exchange(int nlocal, double *buf)
{
  int m,iper,nper,ni;

  k_vectors.sync<LMPHostType>();

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*nlocal;
    for (iper = 0; iper < nper; iper++) h_vectors(m,ni++) = buf[n++];
  }

  k_vectors.modify<LMPHostType>();

  return n;
}
