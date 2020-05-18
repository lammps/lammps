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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "dihedral_opls_kokkos.h"
#include <cmath>
#include <cstdlib>
#include "atom_kokkos.h"
#include "comm.h"
#include "neighbor_kokkos.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001
#define SMALLER   0.00001

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralOPLSKokkos<DeviceType>::DihedralOPLSKokkos(LAMMPS *lmp) : DihedralOPLS(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = DAT::tdual_int_scalar("Dihedral:warning_flag");
  d_warning_flag = k_warning_flag.view<DeviceType>();
  h_warning_flag = k_warning_flag.h_view;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralOPLSKokkos<DeviceType>::~DihedralOPLSKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralOPLSKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"dihedral:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"dihedral:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  k_k1.template sync<DeviceType>();
  k_k2.template sync<DeviceType>();
  k_k3.template sync<DeviceType>();
  k_k4.template sync<DeviceType>();

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_dihedrallist.template sync<DeviceType>();
  dihedrallist = neighborKK->k_dihedrallist.view<DeviceType>();
  int ndihedrallist = neighborKK->ndihedrallist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  h_warning_flag() = 0;
  k_warning_flag.template modify<LMPHostType>();
  k_warning_flag.template sync<DeviceType>();

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralOPLSCompute<1,1> >(0,ndihedrallist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralOPLSCompute<0,1> >(0,ndihedrallist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralOPLSCompute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralOPLSCompute<0,0> >(0,ndihedrallist),*this);
    }
  }

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.template sync<LMPHostType>();
  if (h_warning_flag())
    error->warning(FLERR,"Dihedral problem",0);

  if (eflag_global) energy += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void DihedralOPLSKokkos<DeviceType>::operator()(TagDihedralOPLSCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = dihedrallist(n,0);
  const int i2 = dihedrallist(n,1);
  const int i3 = dihedrallist(n,2);
  const int i4 = dihedrallist(n,3);
  const int type = dihedrallist(n,4);

  // 1st bond

  const F_FLOAT vb1x = x(i1,0) - x(i2,0);
  const F_FLOAT vb1y = x(i1,1) - x(i2,1);
  const F_FLOAT vb1z = x(i1,2) - x(i2,2);

  // 2nd bond

  const F_FLOAT vb2x = x(i3,0) - x(i2,0);
  const F_FLOAT vb2y = x(i3,1) - x(i2,1);
  const F_FLOAT vb2z = x(i3,2) - x(i2,2);

  const F_FLOAT vb2xm = -vb2x;
  const F_FLOAT vb2ym = -vb2y;
  const F_FLOAT vb2zm = -vb2z;

  // 3rd bond

  const F_FLOAT vb3x = x(i4,0) - x(i3,0);
  const F_FLOAT vb3y = x(i4,1) - x(i3,1);
  const F_FLOAT vb3z = x(i4,2) - x(i3,2);

  // c0 calculation

  const F_FLOAT sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  const F_FLOAT sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
  const F_FLOAT sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

  const F_FLOAT rb1 = sqrt(sb1);
  const F_FLOAT rb3 = sqrt(sb3);

  const F_FLOAT c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

  // 1st and 2nd angle

  const F_FLOAT b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
  const F_FLOAT b1mag = sqrt(b1mag2);
  const F_FLOAT b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
  const F_FLOAT b2mag = sqrt(b2mag2);
  const F_FLOAT b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
  const F_FLOAT b3mag = sqrt(b3mag2);

  F_FLOAT ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
  const F_FLOAT r12c1 = 1.0 / (b1mag*b2mag);
  const F_FLOAT c1mag = ctmp * r12c1;

  ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
  const F_FLOAT r12c2 = 1.0 / (b2mag*b3mag);
  const F_FLOAT c2mag = ctmp * r12c2;

  // cos and sin of 2 angles and final c

  F_FLOAT sin2 = MAX(1.0 - c1mag*c1mag,0.0);
  F_FLOAT sc1 = sqrt(sin2);
  if (sc1 < SMALL) sc1 = SMALL;
  sc1 = 1.0/sc1;

  sin2 = MAX(1.0 - c2mag*c2mag,0.0);
  F_FLOAT sc2 = sqrt(sin2);
  if (sc2 < SMALL) sc2 = SMALL;
  sc2 = 1.0/sc2;

  const F_FLOAT s1 = sc1 * sc1;
  const F_FLOAT s2 = sc2 * sc2;
  F_FLOAT s12 = sc1 * sc2;
  F_FLOAT c = (c0 + c1mag*c2mag) * s12;

  const F_FLOAT cx = vb1y*vb2z - vb1z*vb2y;
  const F_FLOAT cy = vb1z*vb2x - vb1x*vb2z;
  const F_FLOAT cz = vb1x*vb2y - vb1y*vb2x;
  const F_FLOAT cmag = sqrt(cx*cx + cy*cy + cz*cz);
  const F_FLOAT dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag;

  // error check

  if ((c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) && !d_warning_flag())
    Kokkos::atomic_fetch_add(&d_warning_flag(),1);

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  // force & energy
  // p = sum (i=1,4) k_i * (1 + (-1)**(i+1)*cos(i*phi) )
  // pd = dp/dc

  F_FLOAT phi = acos(c);
  if (dx < 0.0) phi *= -1.0;
  F_FLOAT si = sin(phi);
  if (fabs(si) < SMALLER) si = SMALLER;
  const F_FLOAT siinv = 1.0/si;

  const F_FLOAT p = d_k1[type]*(1.0 + c) + d_k2[type]*(1.0 - cos(2.0*phi)) +
    d_k3[type]*(1.0 + cos(3.0*phi)) + d_k4[type]*(1.0 - cos(4.0*phi)) ;
  const F_FLOAT pd = d_k1[type] - 2.0*d_k2[type]*sin(2.0*phi)*siinv +
    3.0*d_k3[type]*sin(3.0*phi)*siinv - 4.0*d_k4[type]*sin(4.0*phi)*siinv;

  E_FLOAT edihedral = 0.0;
  if (eflag) edihedral = p;

  const F_FLOAT a = pd;
  c = c * a;
  s12 = s12 * a;
  const F_FLOAT a11 = c*sb1*s1;
  const F_FLOAT a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2));
  const F_FLOAT a33 = c*sb3*s2;
  const F_FLOAT a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12);
  const F_FLOAT a13 = -rb1*rb3*s12;
  const F_FLOAT a23 = r12c2 * (c2mag*c*s2 + c1mag*s12);

  const F_FLOAT sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
  const F_FLOAT sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
  const F_FLOAT sz2  = a12*vb1z + a22*vb2z + a23*vb3z;

  F_FLOAT f1[3],f2[3],f3[3],f4[3];
  f1[0] = a11*vb1x + a12*vb2x + a13*vb3x;
  f1[1] = a11*vb1y + a12*vb2y + a13*vb3y;
  f1[2] = a11*vb1z + a12*vb2z + a13*vb3z;

  f2[0] = -sx2 - f1[0];
  f2[1] = -sy2 - f1[1];
  f2[2] = -sz2 - f1[2];

  f4[0] = a13*vb1x + a23*vb2x + a33*vb3x;
  f4[1] = a13*vb1y + a23*vb2y + a33*vb3y;
  f4[2] = a13*vb1z + a23*vb2z + a33*vb3z;

  f3[0] = sx2 - f4[0];
  f3[1] = sy2 - f4[1];
  f3[2] = sz2 - f4[2];

  // apply force to each of 4 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += f1[0];
    a_f(i1,1) += f1[1];
    a_f(i1,2) += f1[2];
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) += f2[0];
    a_f(i2,1) += f2[1];
    a_f(i2,2) += f2[2];
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += f3[0];
    a_f(i3,1) += f3[1];
    a_f(i3,2) += f3[2];
  }

  if (NEWTON_BOND || i4 < nlocal) {
    a_f(i4,0) += f4[0];
    a_f(i4,1) += f4[1];
    a_f(i4,2) += f4[2];
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,edihedral,f1,f3,f4,
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void DihedralOPLSKokkos<DeviceType>::operator()(TagDihedralOPLSCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralOPLSCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralOPLSKokkos<DeviceType>::allocate()
{
  DihedralOPLS::allocate();

  int n = atom->ndihedraltypes;
  k_k1 = DAT::tdual_ffloat_1d("DihedralOPLS::k1",n+1);
  k_k2 = DAT::tdual_ffloat_1d("DihedralOPLS::k2",n+1);
  k_k3 = DAT::tdual_ffloat_1d("DihedralOPLS::k3",n+1);
  k_k4 = DAT::tdual_ffloat_1d("DihedralOPLS::k4",n+1);

  d_k1 = k_k1.template view<DeviceType>();
  d_k2 = k_k2.template view<DeviceType>();
  d_k3 = k_k3.template view<DeviceType>();
  d_k4 = k_k4.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralOPLSKokkos<DeviceType>::coeff(int narg, char **arg)
{
  DihedralOPLS::coeff(narg, arg);

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_k1.h_view[i] = k1[i];
    k_k2.h_view[i] = k2[i];
    k_k3.h_view[i] = k3[i];
    k_k4.h_view[i] = k4[i];
  }

  k_k1.template modify<LMPHostType>();
  k_k2.template modify<LMPHostType>();
  k_k3.template modify<LMPHostType>();
  k_k4.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralOPLSKokkos<DeviceType>::read_restart(FILE *fp)
{
  DihedralOPLS::read_restart(fp);

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_k1.h_view[i] = k1[i];
    k_k2.h_view[i] = k2[i];
    k_k3.h_view[i] = k3[i];
    k_k4.h_view[i] = k4[i];
  }

  k_k1.template modify<LMPHostType>();
  k_k2.template modify<LMPHostType>();
  k_k3.template modify<LMPHostType>();
  k_k4.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
          = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
KOKKOS_INLINE_FUNCTION
void DihedralOPLSKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                        F_FLOAT &edihedral, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                        const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                        const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                        const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const
{
  E_FLOAT edihedralquarter;
  F_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.view<DeviceType>();

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += edihedral;
      else {
        edihedralquarter = 0.25*edihedral;
        if (i1 < nlocal) ev.evdwl += edihedralquarter;
        if (i2 < nlocal) ev.evdwl += edihedralquarter;
        if (i3 < nlocal) ev.evdwl += edihedralquarter;
        if (i4 < nlocal) ev.evdwl += edihedralquarter;
      }
    }
    if (eflag_atom) {
      edihedralquarter = 0.25*edihedral;
      if (newton_bond || i1 < nlocal) v_eatom[i1] += edihedralquarter;
      if (newton_bond || i2 < nlocal) v_eatom[i2] += edihedralquarter;
      if (newton_bond || i3 < nlocal) v_eatom[i3] += edihedralquarter;
      if (newton_bond || i4 < nlocal) v_eatom[i4] += edihedralquarter;
    }
  }

  if (vflag_either) {
    v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0];
    v[1] = vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1];
    v[2] = vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2];
    v[3] = vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1];
    v[4] = vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2];
    v[5] = vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2];

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i1 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
        if (i2 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
        if (i3 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
        if (i4 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i1 < nlocal) {
        v_vatom(i1,0) += 0.25*v[0];
        v_vatom(i1,1) += 0.25*v[1];
        v_vatom(i1,2) += 0.25*v[2];
        v_vatom(i1,3) += 0.25*v[3];
        v_vatom(i1,4) += 0.25*v[4];
        v_vatom(i1,5) += 0.25*v[5];
      }
      if (newton_bond || i2 < nlocal) {
        v_vatom(i2,0) += 0.25*v[0];
        v_vatom(i2,1) += 0.25*v[1];
        v_vatom(i2,2) += 0.25*v[2];
        v_vatom(i2,3) += 0.25*v[3];
        v_vatom(i2,4) += 0.25*v[4];
        v_vatom(i2,5) += 0.25*v[5];
      }
      if (newton_bond || i3 < nlocal) {
        v_vatom(i3,0) += 0.25*v[0];
        v_vatom(i3,1) += 0.25*v[1];
        v_vatom(i3,2) += 0.25*v[2];
        v_vatom(i3,3) += 0.25*v[3];
        v_vatom(i3,4) += 0.25*v[4];
        v_vatom(i3,5) += 0.25*v[5];
      }
      if (newton_bond || i4 < nlocal) {
        v_vatom(i4,0) += 0.25*v[0];
        v_vatom(i4,1) += 0.25*v[1];
        v_vatom(i4,2) += 0.25*v[2];
        v_vatom(i4,3) += 0.25*v[3];
        v_vatom(i4,4) += 0.25*v[4];
        v_vatom(i4,5) += 0.25*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class DihedralOPLSKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class DihedralOPLSKokkos<LMPHostType>;
#endif
}

