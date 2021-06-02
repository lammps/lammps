// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (Materials Design)
------------------------------------------------------------------------- */

#include "improper_class2_kokkos.h"
#include <cmath>
#include "atom_kokkos.h"
#include "neighbor_kokkos.h"
#include "force.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperClass2Kokkos<DeviceType>::ImproperClass2Kokkos(LAMMPS *lmp) : ImproperClass2(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = DAT::tdual_int_scalar("Dihedral:warning_flag");
  d_warning_flag = k_warning_flag.view<DeviceType>();
  h_warning_flag = k_warning_flag.h_view;

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperClass2Kokkos<DeviceType>::~ImproperClass2Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperClass2Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    //if(k_eatom.extent(0)<maxeatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"improper:eatom");
      d_eatom = k_eatom.template view<DeviceType>();
    //}
  }
  if (vflag_atom) {
    //if(k_vatom.extent(0)<maxvatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"improper:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    //}
  }

  //atomKK->sync(execution_space,datamask_read);
  k_k0.template sync<DeviceType>();
  k_chi0.template sync<DeviceType>();
  k_aa_k1.template sync<DeviceType>();
  k_aa_k2.template sync<DeviceType>();
  k_aa_k3.template sync<DeviceType>();
  k_aa_theta0_1.template sync<DeviceType>();
  k_aa_theta0_2.template sync<DeviceType>();
  k_aa_theta0_3 .template sync<DeviceType>();
  k_setflag.template sync<DeviceType>();
  k_setflag_i.template sync<DeviceType>();
  k_setflag_aa.template sync<DeviceType>();

  //if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  //else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_improperlist.template sync<DeviceType>();
  improperlist = neighborKK->k_improperlist.view<DeviceType>();
  int nimproperlist = neighborKK->nimproperlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  h_warning_flag() = 0;
  k_warning_flag.template modify<LMPHostType>();
  k_warning_flag.template sync<DeviceType>();

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  // Improper energy/force

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperClass2Compute<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperClass2Compute<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperClass2Compute<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperClass2Compute<0,0> >(0,nimproperlist),*this);
    }
  }
  if (eflag_global) energy += ev.evdwl;

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.template sync<LMPHostType>();
  if (h_warning_flag())
    error->warning(FLERR,"Improper problem");

  // Angle-Angle energy/force

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperClass2AngleAngle<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperClass2AngleAngle<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperClass2AngleAngle<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperClass2AngleAngle<0,0> >(0,nimproperlist),*this);
    }
  }
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void ImproperClass2Kokkos<DeviceType>::operator()(TagImproperClass2Compute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  int i, j, k;
  F_FLOAT delr[3][3],rmag[3],rinvmag[3],rmag2[3];
  F_FLOAT theta[3],costheta[3],sintheta[3];
  F_FLOAT cossqtheta[3],sinsqtheta[3],invstheta[3];
  F_FLOAT rABxrCB[3],rDBxrAB[3],rCBxrDB[3];
  F_FLOAT ddelr[3][4],dr[3][4][3],dinvr[3][4][3];
  F_FLOAT dthetadr[3][4][3],dinvsth[3][4][3];
  F_FLOAT dinv3r[4][3],dinvs3r[3][4][3];
  F_FLOAT drCBxrDB[3],rCBxdrDB[3],drDBxrAB[3],rDBxdrAB[3];
  F_FLOAT drABxrCB[3],rABxdrCB[3];
  F_FLOAT dot1,dot2,dd[3];
  F_FLOAT fdot[3][4][3],ftmp,invs3r[3],inv3r;
  F_FLOAT drAB[3][4][3],drCB[3][4][3],drDB[3][4][3];
  F_FLOAT dchi[3][4][3],dtotalchi[4][3];
  F_FLOAT fabcd[4][3];

  F_FLOAT t,tt1,tt3,sc1;
  F_FLOAT dotCBDBAB,dotDBABCB,dotABCBDB;
  F_FLOAT schiABCD,chiABCD,schiCBDA,chiCBDA,schiDBAC,chiDBAC;
  F_FLOAT chi,deltachi,d2chi,cossin2;
  F_FLOAT eimproper;

  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  if (d_k0[type] != 0.0) {

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
        for (k = 0; k < 3; k++) {
          dthetadr[i][j][k] = 0.0;
          drAB[i][j][k] = 0.0;
          drCB[i][j][k] = 0.0;
          drDB[i][j][k] = 0.0;
        }

    // difference vectors

    delr[0][0] = x(i1,0) - x(i2,0);
    delr[0][1] = x(i1,1) - x(i2,1);
    delr[0][2] = x(i1,2) - x(i2,2);

    delr[1][0] = x(i3,0) - x(i2,0);
    delr[1][1] = x(i3,1) - x(i2,1);
    delr[1][2] = x(i3,2) - x(i2,2);

    delr[2][0] = x(i4,0) - x(i2,0);
    delr[2][1] = x(i4,1) - x(i2,1);
    delr[2][2] = x(i4,2) - x(i2,2);

    // bond lengths and associated values

    for (i = 0; i < 3; i++) {
      rmag2[i] = delr[i][0]*delr[i][0] + delr[i][1]*delr[i][1] + delr[i][2]*delr[i][2];
      rmag[i] = sqrt(rmag2[i]);
      rinvmag[i] = 1.0/rmag[i];
    }

    // angle ABC, CBD, ABD

    costheta[0] = (delr[0][0]*delr[1][0] + delr[0][1]*delr[1][1] +
                   delr[0][2]*delr[1][2]) / (rmag[0]*rmag[1]);
    costheta[1] = (delr[1][0]*delr[2][0] + delr[1][1]*delr[2][1] +
                   delr[1][2]*delr[2][2]) / (rmag[1]*rmag[2]);
    costheta[2] = (delr[0][0]*delr[2][0] + delr[0][1]*delr[2][1] +
                   delr[0][2]*delr[2][2]) / (rmag[0]*rmag[2]);

    // sin and cos of improper

    F_FLOAT s1 = 1.0 - costheta[1]*costheta[1];
    if (s1 < SMALL) s1 = SMALL;
    s1 = 1.0 / s1;

    F_FLOAT s2 = 1.0 - costheta[2]*costheta[2];
    if (s2 < SMALL) s2 = SMALL;
    s2 = 1.0 / s2;

    F_FLOAT s12 = sqrt(s1*s2);
    F_FLOAT c = (costheta[1]*costheta[2] + costheta[0]) * s12;

    // error check

    /*
    if ((c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) && !d_warning_flag())
      Kokkos::atomic_fetch_add(&d_warning_flag(),1);
    */
    if ((costheta[0] == -1.0 || costheta[1] == -1.0 || costheta[2] == -1.0) && !d_warning_flag())
      Kokkos::atomic_fetch_add(&d_warning_flag(),1);

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    F_FLOAT s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;

    for (i = 0; i < 3; i++) {
      if (costheta[i] > 1.0)  costheta[i] = 1.0;
      if (costheta[i] < -1.0) costheta[i] = -1.0;
      theta[i] = acos(costheta[i]);
      cossqtheta[i] = costheta[i]*costheta[i];
      sintheta[i] = sin(theta[i]);
      invstheta[i] = 1.0/sintheta[i];
      sinsqtheta[i] = sintheta[i]*sintheta[i];
    }

    // cross & dot products

    //cross(delr[0],delr[1],rABxrCB);
    rABxrCB[0] = delr[0][1]*delr[1][2] - delr[0][2]*delr[1][1];
    rABxrCB[1] = delr[0][2]*delr[1][0] - delr[0][0]*delr[1][2];
    rABxrCB[2] = delr[0][0]*delr[1][1] - delr[0][1]*delr[1][0];

    //cross(delr[2],delr[0],rDBxrAB);
    rDBxrAB[0] = delr[2][1]*delr[0][2] - delr[2][2]*delr[0][1];
    rDBxrAB[1] = delr[2][2]*delr[0][0] - delr[2][0]*delr[0][2];
    rDBxrAB[2] = delr[2][0]*delr[0][1] - delr[2][1]*delr[0][0];

    //cross(delr[1],delr[2],rCBxrDB);
    rCBxrDB[0] = delr[1][1]*delr[2][2] - delr[1][2]*delr[2][1];
    rCBxrDB[1] = delr[1][2]*delr[2][0] - delr[1][0]*delr[2][2];
    rCBxrDB[2] = delr[1][0]*delr[2][1] - delr[1][1]*delr[2][0];

    //dotCBDBAB = dot(rCBxrDB,delr[0]);
    dotCBDBAB = rCBxrDB[0]*delr[0][0] + rCBxrDB[1]*delr[0][1] + rCBxrDB[2]*delr[0][2];

    //dotDBABCB = dot(rDBxrAB,delr[1]);
    dotDBABCB = rDBxrAB[0]*delr[1][0] + rDBxrAB[1]*delr[1][1] + rDBxrAB[2]*delr[1][2];

    //dotABCBDB = dot(rABxrCB,delr[2]);
    dotABCBDB = rABxrCB[0]*delr[2][0] + rABxrCB[1]*delr[2][1] + rABxrCB[2]*delr[2][2];

    t = rmag[0] * rmag[1] * rmag[2];
    inv3r = 1.0/t;
    invs3r[0] = invstheta[1] * inv3r;
    invs3r[1] = invstheta[2] * inv3r;
    invs3r[2] = invstheta[0] * inv3r;

    // chi ABCD, CBDA, DBAC: final chi is average of three

    schiABCD = dotCBDBAB * invs3r[0];
    chiABCD = asin(schiABCD);
    schiCBDA = dotDBABCB * invs3r[1];
    chiCBDA = asin(schiCBDA);
    schiDBAC = dotABCBDB * invs3r[2];
    chiDBAC = asin(schiDBAC);

    chi = (chiABCD + chiCBDA + chiDBAC) / 3.0;
    deltachi = chi - d_chi0[type];
    d2chi = deltachi * deltachi;

    // energy

    if (eflag) eimproper = d_k0[type]*d2chi;

    // forces

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
        ddelr[i][j] = 0.0;

    ddelr[0][0] = 1.0;
    ddelr[0][1] = -1.0;
    ddelr[1][1] = -1.0;
    ddelr[1][2] = 1.0;
    ddelr[2][1] = -1.0;
    ddelr[2][3] = 1.0;

    // compute d(|r|)/dr and d(1/|r|)/dr for each direction, bond and atom

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
        for (k = 0; k < 3; k++) {
          dr[i][j][k] = delr[i][k] * ddelr[i][j] / rmag[i];
          dinvr[i][j][k] = -dr[i][j][k] / rmag2[i];
        }

    // compute d(1 / (|r_AB| * |r_CB| * |r_DB|) / dr

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        dinv3r[i][j] = rinvmag[1] * (rinvmag[2] * dinvr[0][i][j] +
                                     rinvmag[0] * dinvr[2][i][j]) +
          rinvmag[2] * rinvmag[0] * dinvr[1][i][j];

    // compute d(theta)/d(r) for 3 angles
    // angleABC

    tt1 = costheta[0] / rmag2[0];
    tt3 = costheta[0] / rmag2[1];
    sc1 = 1.0 / sqrt(1.0 - cossqtheta[0]);

    dthetadr[0][0][0] = sc1 * ((tt1 * delr[0][0]) -
                               (delr[1][0] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][0][1] = sc1 * ((tt1 * delr[0][1]) -
                               (delr[1][1] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][0][2] = sc1 * ((tt1 * delr[0][2]) -
                               (delr[1][2] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][1][0] = -sc1 * ((tt1 * delr[0][0]) -
                                (delr[1][0] * rinvmag[0] * rinvmag[1]) +
                                (tt3 * delr[1][0]) -
                                (delr[0][0] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][1][1] = -sc1 * ((tt1 * delr[0][1]) -
                                (delr[1][1] * rinvmag[0] * rinvmag[1]) +
                                (tt3 * delr[1][1]) -
                                (delr[0][1] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][1][2] = -sc1 * ((tt1 * delr[0][2]) -
                                (delr[1][2] * rinvmag[0] * rinvmag[1]) +
                                (tt3 * delr[1][2]) -
                                (delr[0][2] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][2][0] = sc1 * ((tt3 * delr[1][0]) -
                               (delr[0][0] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][2][1] = sc1 * ((tt3 * delr[1][1]) -
                               (delr[0][1] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][2][2] = sc1 * ((tt3 * delr[1][2]) -
                               (delr[0][2] * rinvmag[0] * rinvmag[1]));


    // angleCBD

    tt1 = costheta[1] / rmag2[1];
    tt3 = costheta[1] / rmag2[2];
    sc1 = 1.0 / sqrt(1.0 - cossqtheta[1]);

    dthetadr[1][2][0] = sc1 * ((tt1 * delr[1][0]) -
                               (delr[2][0] * rinvmag[1] * rinvmag[2]));
    dthetadr[1][2][1] = sc1 * ((tt1 * delr[1][1]) -
                               (delr[2][1] * rinvmag[1] * rinvmag[2]));
    dthetadr[1][2][2] = sc1 * ((tt1 * delr[1][2]) -
                               (delr[2][2] * rinvmag[1] * rinvmag[2]));
    dthetadr[1][1][0] = -sc1 * ((tt1 * delr[1][0]) -
                                (delr[2][0] * rinvmag[1] * rinvmag[2]) +
                                (tt3 * delr[2][0]) -
                                (delr[1][0] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][1][1] = -sc1 * ((tt1 * delr[1][1]) -
                                (delr[2][1] * rinvmag[1] * rinvmag[2]) +
                                (tt3 * delr[2][1]) -
                                (delr[1][1] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][1][2] = -sc1 * ((tt1 * delr[1][2]) -
                                (delr[2][2] * rinvmag[1] * rinvmag[2]) +
                                (tt3 * delr[2][2]) -
                                (delr[1][2] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][3][0] = sc1 * ((tt3 * delr[2][0]) -
                               (delr[1][0] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][3][1] = sc1 * ((tt3 * delr[2][1]) -
                               (delr[1][1] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][3][2] = sc1 * ((tt3 * delr[2][2]) -
                               (delr[1][2] * rinvmag[2] * rinvmag[1]));

    // angleABD

    tt1 = costheta[2] / rmag2[0];
    tt3 = costheta[2] / rmag2[2];
    sc1 = 1.0 / sqrt(1.0 - cossqtheta[2]);

    dthetadr[2][0][0] = sc1 * ((tt1 * delr[0][0]) -
                               (delr[2][0] * rinvmag[0] * rinvmag[2]));
    dthetadr[2][0][1] = sc1 * ((tt1 * delr[0][1]) -
                               (delr[2][1] * rinvmag[0] * rinvmag[2]));
    dthetadr[2][0][2] = sc1 * ((tt1 * delr[0][2]) -
                               (delr[2][2] * rinvmag[0] * rinvmag[2]));
    dthetadr[2][1][0] = -sc1 * ((tt1 * delr[0][0]) -
                                (delr[2][0] * rinvmag[0] * rinvmag[2]) +
                                (tt3 * delr[2][0]) -
                                (delr[0][0] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][1][1] = -sc1 * ((tt1 * delr[0][1]) -
                                (delr[2][1] * rinvmag[0] * rinvmag[2]) +
                                (tt3 * delr[2][1]) -
                                (delr[0][1] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][1][2] = -sc1 * ((tt1 * delr[0][2]) -
                                (delr[2][2] * rinvmag[0] * rinvmag[2]) +
                                (tt3 * delr[2][2]) -
                                (delr[0][2] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][3][0] = sc1 * ((tt3 * delr[2][0]) -
                               (delr[0][0] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][3][1] = sc1 * ((tt3 * delr[2][1]) -
                               (delr[0][1] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][3][2] = sc1 * ((tt3 * delr[2][2]) -
                               (delr[0][2] * rinvmag[2] * rinvmag[0]));

    // compute d( 1 / sin(theta))/dr
    // i = angle, j = atom, k = direction

    for (i = 0; i < 3; i++) {
      cossin2 = -costheta[i] / sinsqtheta[i];
      for (j = 0; j < 4; j++)
        for (k = 0; k < 3; k++)
          dinvsth[i][j][k] = cossin2 * dthetadr[i][j][k];
    }

    // compute d(1 / sin(theta) * |r_AB| * |r_CB| * |r_DB|)/dr
    // i = angle, j = atom

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++) {
        dinvs3r[0][i][j] = (invstheta[1] * dinv3r[i][j]) +
          (inv3r * dinvsth[1][i][j]);
        dinvs3r[1][i][j] = (invstheta[2] * dinv3r[i][j]) +
          (inv3r * dinvsth[2][i][j]);
        dinvs3r[2][i][j] = (invstheta[0] * dinv3r[i][j]) +
          (inv3r * dinvsth[0][i][j]);
      }

    // drCB(i,j,k), etc
    // i = vector X'/Y'/Z', j = atom A/B/C/D, k = direction X/Y/Z

    for (i = 0; i < 3; i++) {
      drCB[i][1][i] = -1.0;
      drAB[i][1][i] = -1.0;
      drDB[i][1][i] = -1.0;
      drDB[i][3][i] = 1.0;
      drCB[i][2][i] = 1.0;
      drAB[i][0][i] = 1.0;
    }

    // d((r_CB x r_DB) dot r_AB)

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++) {
        //cross(delr[1],drDB[i][j],rCBxdrDB);
        rCBxdrDB[0] = delr[1][1]*drDB[i][j][2] - delr[1][2]*drDB[i][j][1];
        rCBxdrDB[1] = delr[1][2]*drDB[i][j][0] - delr[1][0]*drDB[i][j][2];
        rCBxdrDB[2] = delr[1][0]*drDB[i][j][1] - delr[1][1]*drDB[i][j][0];

        //cross(drCB[i][j],delr[2],drCBxrDB);
        drCBxrDB[0] = drCB[i][j][1]*delr[2][2] - drCB[i][j][2]*delr[2][1];
        drCBxrDB[1] = drCB[i][j][2]*delr[2][0] - drCB[i][j][0]*delr[2][2];
        drCBxrDB[2] = drCB[i][j][0]*delr[2][1] - drCB[i][j][1]*delr[2][0];

        for (k = 0; k < 3; k++) dd[k] = rCBxdrDB[k] + drCBxrDB[k];
        //dot1 = dot(dd,delr[0]);
        dot1 = dd[0]*delr[0][0] + dd[1]*delr[0][1] + dd[2]*delr[0][2];

        //dot2 = dot(rCBxrDB,drAB[i][j]);
        dot2 = rCBxrDB[0]*drAB[i][j][0] + rCBxrDB[1]*drAB[i][j][1] + rCBxrDB[2]*drAB[i][j][2];

        fdot[0][j][i] = dot1 + dot2;
      }

    // d((r_DB x r_AB) dot r_CB)

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++) {
        //cross(delr[2],drAB[i][j],rDBxdrAB);
        rDBxdrAB[0] = delr[2][1]*drAB[i][j][2] - delr[2][2]*drAB[i][j][1];
        rDBxdrAB[1] = delr[2][2]*drAB[i][j][0] - delr[2][0]*drAB[i][j][2];
        rDBxdrAB[2] = delr[2][0]*drAB[i][j][1] - delr[2][1]*drAB[i][j][0];

        //cross(drDB[i][j],delr[0],drDBxrAB);
        drDBxrAB[0] = drDB[i][j][1]*delr[0][2] - drDB[i][j][2]*delr[0][1];
        drDBxrAB[1] = drDB[i][j][2]*delr[0][0] - drDB[i][j][0]*delr[0][2];
        drDBxrAB[2] = drDB[i][j][0]*delr[0][1] - drDB[i][j][1]*delr[0][0];

        for (k = 0; k < 3; k++) dd[k] = rDBxdrAB[k] + drDBxrAB[k];

        //dot1 = dot(dd,delr[1]);
        dot1 = dd[0]*delr[1][0] + dd[1]*delr[1][1] + dd[2]*delr[1][2];

        //dot2 = dot(rDBxrAB,drCB[i][j]);
        dot2 = rDBxrAB[0]*drCB[i][j][0] + rDBxrAB[1]*drCB[i][j][1] + rDBxrAB[2]*drCB[i][j][2];

        fdot[1][j][i] = dot1 + dot2;
      }

    // d((r_AB x r_CB) dot r_DB)

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++) {
        //cross(delr[0],drCB[i][j],rABxdrCB);
        rABxdrCB[0] = delr[0][1]*drCB[i][j][2] - delr[0][2]*drCB[i][j][1];
        rABxdrCB[1] = delr[0][2]*drCB[i][j][0] - delr[0][0]*drCB[i][j][2];
        rABxdrCB[2] = delr[0][0]*drCB[i][j][1] - delr[0][1]*drCB[i][j][0];

        //cross(drAB[i][j],delr[1],drABxrCB);
        drABxrCB[0] = drAB[i][j][1]*delr[1][2] - drAB[i][j][2]*delr[1][1];
        drABxrCB[1] = drAB[i][j][2]*delr[1][0] - drAB[i][j][0]*delr[1][2];
        drABxrCB[2] = drAB[i][j][0]*delr[1][1] - drAB[i][j][1]*delr[1][0];

        for (k = 0; k < 3; k++) dd[k] = rABxdrCB[k] + drABxrCB[k];

        //dot1 = dot(dd,delr[2]);
        dot1 = dd[0]*delr[2][0] + dd[1]*delr[2][1] + dd[2]*delr[2][2];

        //dot2 = dot(rABxrCB,drDB[i][j]);
        dot2 = rABxrCB[0]*drDB[i][j][0] +rABxrCB[1]*drDB[i][j][1] +rABxrCB[2]*drDB[i][j][2];

        fdot[2][j][i] = dot1 + dot2;
      }

    // force on each atom

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++) {
        ftmp = (fdot[0][i][j] * invs3r[0]) + (dinvs3r[0][i][j] * dotCBDBAB);
        dchi[0][i][j] = ftmp / cos(chiABCD);
        ftmp = (fdot[1][i][j] * invs3r[1]) + (dinvs3r[1][i][j] * dotDBABCB);
        dchi[1][i][j] = ftmp / cos(chiCBDA);
        ftmp = (fdot[2][i][j] * invs3r[2]) + (dinvs3r[2][i][j] * dotABCBDB);
        dchi[2][i][j] = ftmp / cos(chiDBAC);
        dtotalchi[i][j] = (dchi[0][i][j]+dchi[1][i][j]+dchi[2][i][j]) / 3.0;
      }

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        fabcd[i][j] = -2.0*d_k0[type] * deltachi*dtotalchi[i][j];

    // apply force to each of 4 atoms

    F_FLOAT f1[3],f2[3],f3[3],f4[3];

    for (i = 0; i < 3; i++) {
      f1[i] = fabcd[0][i];
      f2[i] = fabcd[1][i];
      f3[i] = fabcd[2][i];
      f4[i] = fabcd[3][i];
    }

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
      ev_tally(ev,i1,i2,i3,i4,eimproper,f1,f3,f4,
               delr[0][0], delr[0][1], delr[0][2],
               delr[1][0], delr[1][1], delr[1][2],
               delr[2][0]- delr[1][0], delr[2][1]-delr[1][1], delr[2][2]-delr[1][2]);
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void ImproperClass2Kokkos<DeviceType>::operator()(TagImproperClass2Compute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperClass2Compute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void ImproperClass2Kokkos<DeviceType>::operator()(TagImproperClass2AngleAngle<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  int i,j,k;
  F_FLOAT eimproper;
  F_FLOAT delxAB,delyAB,delzAB,rABmag2,rAB;
  F_FLOAT delxBC,delyBC,delzBC,rBCmag2,rBC;
  F_FLOAT delxBD,delyBD,delzBD,rBDmag2,rBD;
  F_FLOAT costhABC,thetaABC,costhABD;
  F_FLOAT thetaABD,costhCBD,thetaCBD,dthABC,dthCBD,dthABD;
  F_FLOAT sc1,t1,t3,r12;
  F_FLOAT dthetadr[3][4][3],fabcd[4][3];

  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  if ((d_aa_k1[type] != 0.0) || (d_aa_k2[type] != 0.0) || (d_aa_k3[type] != 0.0)) {

    // difference vectors

    delxAB = x(i1,0) - x(i2,0);
    delyAB = x(i1,1) - x(i2,1);
    delzAB = x(i1,2) - x(i2,2);

    delxBC = x(i3,0) - x(i2,0);
    delyBC = x(i3,1) - x(i2,1);
    delzBC = x(i3,2) - x(i2,2);

    delxBD = x(i4,0) - x(i2,0);
    delyBD = x(i4,1) - x(i2,1);
    delzBD = x(i4,2) - x(i2,2);

    // bond lengths

    rABmag2 = delxAB*delxAB + delyAB*delyAB + delzAB*delzAB;
    rAB = sqrt(rABmag2);
    rBCmag2 = delxBC*delxBC + delyBC*delyBC + delzBC*delzBC;
    rBC = sqrt(rBCmag2);
    rBDmag2 = delxBD*delxBD + delyBD*delyBD + delzBD*delzBD;
    rBD = sqrt(rBDmag2);

    // angle ABC, ABD, CBD

    costhABC = (delxAB*delxBC + delyAB*delyBC + delzAB*delzBC) / (rAB * rBC);
    if (costhABC > 1.0)  costhABC = 1.0;
    if (costhABC < -1.0) costhABC = -1.0;
    thetaABC = acos(costhABC);

    costhABD = (delxAB*delxBD + delyAB*delyBD + delzAB*delzBD) / (rAB * rBD);
    if (costhABD > 1.0)  costhABD = 1.0;
    if (costhABD < -1.0) costhABD = -1.0;
    thetaABD = acos(costhABD);

    costhCBD = (delxBC*delxBD + delyBC*delyBD + delzBC*delzBD) /(rBC * rBD);
    if (costhCBD > 1.0)  costhCBD = 1.0;
    if (costhCBD < -1.0) costhCBD = -1.0;
    thetaCBD = acos(costhCBD);

    dthABC = thetaABC - d_aa_theta0_1[type];
    dthABD = thetaABD - d_aa_theta0_2[type];
    dthCBD = thetaCBD - d_aa_theta0_3[type];

    // energy

    if (eflag) eimproper = d_aa_k2[type] * dthABC * dthABD +
                 d_aa_k1[type] * dthABC * dthCBD +
                 d_aa_k3[type] * dthABD * dthCBD;

    // d(theta)/d(r) array

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
        for (k = 0; k < 3; k++)
          dthetadr[i][j][k] = 0.0;

    // angle ABC

    sc1 = sqrt(1.0/(1.0 - costhABC*costhABC));
    t1 = costhABC / rABmag2;
    t3 = costhABC / rBCmag2;
    r12 = 1.0 / (rAB * rBC);

    dthetadr[0][0][0] = sc1 * ((t1 * delxAB) - (delxBC * r12));
    dthetadr[0][0][1] = sc1 * ((t1 * delyAB) - (delyBC * r12));
    dthetadr[0][0][2] = sc1 * ((t1 * delzAB) - (delzBC * r12));
    dthetadr[0][1][0] = sc1 * ((-t1 * delxAB) + (delxBC * r12) +
                               (-t3 * delxBC) + (delxAB * r12));
    dthetadr[0][1][1] = sc1 * ((-t1 * delyAB) + (delyBC * r12) +
                               (-t3 * delyBC) + (delyAB * r12));
    dthetadr[0][1][2] = sc1 * ((-t1 * delzAB) + (delzBC * r12) +
                               (-t3 * delzBC) + (delzAB * r12));
    dthetadr[0][2][0] = sc1 * ((t3 * delxBC) - (delxAB * r12));
    dthetadr[0][2][1] = sc1 * ((t3 * delyBC) - (delyAB * r12));
    dthetadr[0][2][2] = sc1 * ((t3 * delzBC) - (delzAB * r12));

    // angle CBD

    sc1 = sqrt(1.0/(1.0 - costhCBD*costhCBD));
    t1 = costhCBD / rBCmag2;
    t3 = costhCBD / rBDmag2;
    r12 = 1.0 / (rBC * rBD);

    dthetadr[1][2][0] = sc1 * ((t1 * delxBC) - (delxBD * r12));
    dthetadr[1][2][1] = sc1 * ((t1 * delyBC) - (delyBD * r12));
    dthetadr[1][2][2] = sc1 * ((t1 * delzBC) - (delzBD * r12));
    dthetadr[1][1][0] = sc1 * ((-t1 * delxBC) + (delxBD * r12) +
                               (-t3 * delxBD) + (delxBC * r12));
    dthetadr[1][1][1] = sc1 * ((-t1 * delyBC) + (delyBD * r12) +
                               (-t3 * delyBD) + (delyBC * r12));
    dthetadr[1][1][2] = sc1 * ((-t1 * delzBC) + (delzBD * r12) +
                               (-t3 * delzBD) + (delzBC * r12));
    dthetadr[1][3][0] = sc1 * ((t3 * delxBD) - (delxBC * r12));
    dthetadr[1][3][1] = sc1 * ((t3 * delyBD) - (delyBC * r12));
    dthetadr[1][3][2] = sc1 * ((t3 * delzBD) - (delzBC * r12));

    // angle ABD

    sc1 = sqrt(1.0/(1.0 - costhABD*costhABD));
    t1 = costhABD / rABmag2;
    t3 = costhABD / rBDmag2;
    r12 = 1.0 / (rAB * rBD);

    dthetadr[2][0][0] = sc1 * ((t1 * delxAB) - (delxBD * r12));
    dthetadr[2][0][1] = sc1 * ((t1 * delyAB) - (delyBD * r12));
    dthetadr[2][0][2] = sc1 * ((t1 * delzAB) - (delzBD * r12));
    dthetadr[2][1][0] = sc1 * ((-t1 * delxAB) + (delxBD * r12) +
                               (-t3 * delxBD) + (delxAB * r12));
    dthetadr[2][1][1] = sc1 * ((-t1 * delyAB) + (delyBD * r12) +
                               (-t3 * delyBD) + (delyAB * r12));
    dthetadr[2][1][2] = sc1 * ((-t1 * delzAB) + (delzBD * r12) +
                               (-t3 * delzBD) + (delzAB * r12));
    dthetadr[2][3][0] = sc1 * ((t3 * delxBD) - (delxAB * r12));
    dthetadr[2][3][1] = sc1 * ((t3 * delyBD) - (delyAB * r12));
    dthetadr[2][3][2] = sc1 * ((t3 * delzBD) - (delzAB * r12));

    // angleangle forces

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        fabcd[i][j] = -
          ((d_aa_k1[type] * (dthABC*dthetadr[1][i][j] + dthCBD*dthetadr[0][i][j])) +
           (d_aa_k2[type] * (dthABC*dthetadr[2][i][j] + dthABD*dthetadr[0][i][j])) +
           (d_aa_k3[type] * (dthABD*dthetadr[1][i][j] + dthCBD*dthetadr[2][i][j])));

    // apply force to each of 4 atoms

    F_FLOAT f1[3],f2[3],f3[3],f4[3];

    for (i = 0; i < 3; i++) {
      f1[i] = fabcd[0][i];
      f2[i] = fabcd[1][i];
      f3[i] = fabcd[2][i];
      f4[i] = fabcd[3][i];
    }

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
      ev_tally(ev,i1,i2,i3,i4,eimproper,
               fabcd[0],fabcd[2],fabcd[3],
               delxAB,delyAB,delzAB,delxBC,delyBC,delzBC,
               delxBD-delxBC,delyBD-delyBC,delzBD-delzBC);
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void ImproperClass2Kokkos<DeviceType>::operator()(TagImproperClass2AngleAngle<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperClass2AngleAngle<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperClass2Kokkos<DeviceType>::allocate()
{
  ImproperClass2::allocate();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void ImproperClass2Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  ImproperClass2::coeff(narg, arg);

  int n = atom->nimpropertypes;
  k_k0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::k0",n+1);
  k_chi0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::chi0",n+1);
  k_aa_k1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_k1",n+1);
  k_aa_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_k2",n+1);
  k_aa_k3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_k3",n+1);
  k_aa_theta0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_theta0_1",n+1);
  k_aa_theta0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_theta0_2",n+1);
  k_aa_theta0_3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_theta0_3",n+1);
  k_setflag = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::setflag",n+1);
  k_setflag_i = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::setflag_i",n+1);
  k_setflag_aa = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::setflag_aa",n+1);

  d_k0 = k_k0.template view<DeviceType>();
  d_chi0 = k_chi0.template view<DeviceType>();
  d_aa_k1 = k_aa_k1.template view<DeviceType>();
  d_aa_k2 = k_aa_k2.template view<DeviceType>();
  d_aa_k3 = k_aa_k3.template view<DeviceType>();
  d_aa_theta0_1 = k_aa_theta0_1.template view<DeviceType>();
  d_aa_theta0_2 = k_aa_theta0_2.template view<DeviceType>();
  d_aa_theta0_3 = k_aa_theta0_3.template view<DeviceType>();
  d_setflag = k_setflag.template view<DeviceType>();
  d_setflag_i = k_setflag_i.template view<DeviceType>();
  d_setflag_aa = k_setflag_aa.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k0.h_view[i] = k0[i];
    k_chi0.h_view[i] = chi0[i];
    k_aa_k1.h_view[i] = aa_k1[i];
    k_aa_k2.h_view[i] = aa_k2[i];
    k_aa_k3.h_view[i] = aa_k3[i];
    k_aa_theta0_1.h_view[i] = aa_theta0_1[i];
    k_aa_theta0_2.h_view[i] = aa_theta0_2[i];
    k_aa_theta0_3.h_view[i] = aa_theta0_3[i];
    k_setflag.h_view[i] = setflag[i];
    k_setflag_i.h_view[i] = setflag_i[i];
    k_setflag_aa.h_view[i] = setflag_aa[i];
  }

  k_k0.template modify<LMPHostType>();
  k_chi0.template modify<LMPHostType>();
  k_aa_k1.template modify<LMPHostType>();
  k_aa_k2.template modify<LMPHostType>();
  k_aa_k3.template modify<LMPHostType>();
  k_aa_theta0_1.template modify<LMPHostType>();
  k_aa_theta0_2.template modify<LMPHostType>();
  k_aa_theta0_3 .template modify<LMPHostType>();
  k_setflag.template modify<LMPHostType>();
  k_setflag_i.template modify<LMPHostType>();
  k_setflag_aa.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void ImproperClass2Kokkos<DeviceType>::read_restart(FILE *fp)
{
  ImproperClass2::read_restart(fp);

  int n = atom->nimpropertypes;
  k_k0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::k0",n+1);
  k_chi0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::chi0",n+1);
  k_aa_k1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_k1",n+1);
  k_aa_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_k2",n+1);
  k_aa_k3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_k3",n+1);
  k_aa_theta0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_theta0_1",n+1);
  k_aa_theta0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_theta0_2",n+1);
  k_aa_theta0_3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::aa_theta0_3",n+1);
  k_setflag = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::setflag",n+1);
  k_setflag_i = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::setflag_i",n+1);
  k_setflag_aa = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("ImproperClass2::setflag_aa",n+1);

  d_k0 = k_k0.template view<DeviceType>();
  d_chi0 = k_chi0.template view<DeviceType>();
  d_aa_k1 = k_aa_k1.template view<DeviceType>();
  d_aa_k2 = k_aa_k2.template view<DeviceType>();
  d_aa_k3 = k_aa_k3.template view<DeviceType>();
  d_aa_theta0_1 = k_aa_theta0_1.template view<DeviceType>();
  d_aa_theta0_2 = k_aa_theta0_2.template view<DeviceType>();
  d_aa_theta0_3 = k_aa_theta0_3.template view<DeviceType>();
  d_setflag = k_setflag.template view<DeviceType>();
  d_setflag_i = k_setflag_i.template view<DeviceType>();
  d_setflag_aa = k_setflag_aa.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k0.h_view[i] = k0[i];
    k_chi0.h_view[i] = chi0[i];
    k_aa_k1.h_view[i] = aa_k1[i];
    k_aa_k2.h_view[i] = aa_k2[i];
    k_aa_k3.h_view[i] = aa_k3[i];
    k_aa_theta0_1.h_view[i] = aa_theta0_1[i];
    k_aa_theta0_2.h_view[i] = aa_theta0_2[i];
    k_aa_theta0_3.h_view[i] = aa_theta0_3[i];
    k_setflag.h_view[i] = setflag[i];
    k_setflag_i.h_view[i] = setflag_i[i];
    k_setflag_aa.h_view[i] = setflag_aa[i];
  }

  k_k0.template modify<LMPHostType>();
  k_chi0.template modify<LMPHostType>();
  k_aa_k1.template modify<LMPHostType>();
  k_aa_k2.template modify<LMPHostType>();
  k_aa_k3.template modify<LMPHostType>();
  k_aa_theta0_1.template modify<LMPHostType>();
  k_aa_theta0_2.template modify<LMPHostType>();
  k_aa_theta0_3 .template modify<LMPHostType>();
  k_setflag.template modify<LMPHostType>();
  k_setflag_i.template modify<LMPHostType>();
  k_setflag_aa.template modify<LMPHostType>();
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
void ImproperClass2Kokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                        F_FLOAT &eimproper, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                        const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                        const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                        const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const
{
  E_FLOAT eimproperquarter;
  F_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.view<DeviceType>();

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += eimproper;
      else {
        eimproperquarter = 0.25*eimproper;
        if (i1 < nlocal) ev.evdwl += eimproperquarter;
        if (i2 < nlocal) ev.evdwl += eimproperquarter;
        if (i3 < nlocal) ev.evdwl += eimproperquarter;
        if (i4 < nlocal) ev.evdwl += eimproperquarter;
      }
    }
    if (eflag_atom) {
      eimproperquarter = 0.25*eimproper;
      if (newton_bond || i1 < nlocal) v_eatom[i1] += eimproperquarter;
      if (newton_bond || i2 < nlocal) v_eatom[i2] += eimproperquarter;
      if (newton_bond || i3 < nlocal) v_eatom[i3] += eimproperquarter;
      if (newton_bond || i4 < nlocal) v_eatom[i4] += eimproperquarter;
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
template class ImproperClass2Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ImproperClass2Kokkos<LMPHostType>;
#endif
}

