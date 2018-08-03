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
   Contributing author: Ray Shan (Materials Design)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include "dihedral_class2_kokkos.h"
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
DihedralClass2Kokkos<DeviceType>::DihedralClass2Kokkos(LAMMPS *lmp) : DihedralClass2(lmp)
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
DihedralClass2Kokkos<DeviceType>::~DihedralClass2Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralClass2Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (eflag || vflag) ev_setup(eflag,vflag,0);
  else evflag = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"dihedral:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"dihedral:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  //atomKK->sync(execution_space,datamask_read);
  k_k1.template sync<DeviceType>();
  k_k2.template sync<DeviceType>();
  k_k3.template sync<DeviceType>();
  k_phi1.template sync<DeviceType>();
  k_phi2.template sync<DeviceType>();
  k_phi3.template sync<DeviceType>();
  k_mbt_f1.template sync<DeviceType>();
  k_mbt_f2.template sync<DeviceType>();
  k_mbt_f3.template sync<DeviceType>();
  k_mbt_r0.template sync<DeviceType>();
  k_ebt_f1_1.template sync<DeviceType>();
  k_ebt_f2_1.template sync<DeviceType>();
  k_ebt_f3_1.template sync<DeviceType>();
  k_ebt_r0_1.template sync<DeviceType>();
  k_ebt_f1_2.template sync<DeviceType>();
  k_ebt_f2_2.template sync<DeviceType>();
  k_ebt_f3_2.template sync<DeviceType>();
  k_ebt_r0_2.template sync<DeviceType>();
  k_at_f1_1.template sync<DeviceType>();
  k_at_f2_1.template sync<DeviceType>();
  k_at_f3_1.template sync<DeviceType>();
  k_at_f1_2.template sync<DeviceType>();
  k_at_f2_2.template sync<DeviceType>();
  k_at_f3_2.template sync<DeviceType>();
  k_at_theta0_1.template sync<DeviceType>();
  k_at_theta0_2.template sync<DeviceType>();
  k_aat_k.template sync<DeviceType>();
  k_aat_theta0_1.template sync<DeviceType>();
  k_aat_theta0_2.template sync<DeviceType>();
  k_bb13t_k.template sync<DeviceType>();
  k_bb13t_r10.template sync<DeviceType>();
  k_bb13t_r30.template sync<DeviceType>();
  k_setflag_d.template sync<DeviceType>();
  k_setflag_mbt.template sync<DeviceType>();
  k_setflag_ebt.template sync<DeviceType>();
  k_setflag_at.template sync<DeviceType>();
  k_setflag_aat.template sync<DeviceType>();
  k_setflag_bb13t.template sync<DeviceType>();

  //if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  //else atomKK->modified(execution_space,F_MASK);

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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralClass2Compute<1,1> >(0,ndihedrallist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralClass2Compute<0,1> >(0,ndihedrallist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralClass2Compute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralClass2Compute<0,0> >(0,ndihedrallist),*this);
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void DihedralClass2Kokkos<DeviceType>::operator()(TagDihedralClass2Compute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

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

  // distance: c0 calculation

  const F_FLOAT r1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
  const F_FLOAT r1 = sqrt(r1mag2);
  const F_FLOAT r2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
  const F_FLOAT r2 = sqrt(r2mag2);
  const F_FLOAT r3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
  const F_FLOAT r3 = sqrt(r3mag2);

  const F_FLOAT sb1 = 1.0/r1mag2;
  const F_FLOAT rb1 = 1.0/r1;
  const F_FLOAT sb2 = 1.0/r2mag2;
  const F_FLOAT rb2 = 1.0/r2;
  const F_FLOAT sb3 = 1.0/r3mag2;
  const F_FLOAT rb3 = 1.0/r3;

  const F_FLOAT c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

  // 1st and 2nd angle

  const F_FLOAT r12c1 = rb1*rb2;
  const F_FLOAT r12c2 = rb2*rb3;
  const F_FLOAT costh12 = (vb1x*vb2x + vb1y*vb2y + vb1z*vb2z) * r12c1;
  const F_FLOAT costh13 = c0;
  const F_FLOAT costh23 = (vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z) * r12c2;

  // cos and sin of 2 angles and final c

  F_FLOAT sin2 = MAX(1.0 - costh12*costh12,0.0);
  F_FLOAT sc1 = sqrt(sin2);
  if (sc1 < SMALL) sc1 = SMALL;
  sc1 = 1.0/sc1;

  sin2 = MAX(1.0 - costh23*costh23,0.0);
  F_FLOAT sc2 = sqrt(sin2);
  if (sc2 < SMALL) sc2 = SMALL;
  sc2 = 1.0/sc2;

  const F_FLOAT s1 = sc1 * sc1;
  const F_FLOAT s2 = sc2 * sc2;
  const F_FLOAT s12 = sc1 * sc2;
  F_FLOAT c = (c0 + costh12*costh23) * s12;

  // error check

  if ((c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) && !d_warning_flag())
    Kokkos::atomic_fetch_add(&d_warning_flag(),1);

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  const F_FLOAT cosphi = c;
  F_FLOAT phi = acos(c);

  F_FLOAT sinphi = sqrt(1.0 - c*c);
  sinphi = MAX(sinphi,SMALL);

  // n123 = vb1 x vb2

  const F_FLOAT n123x = vb1y*vb2z - vb1z*vb2y;
  const F_FLOAT n123y = vb1z*vb2x - vb1x*vb2z;
  const F_FLOAT n123z = vb1x*vb2y - vb1y*vb2x;
  const F_FLOAT n123_dot_vb3 = n123x*vb3x + n123y*vb3y + n123z*vb3z;
  if (n123_dot_vb3 > 0.0) {
    phi = -phi;
    sinphi = -sinphi;
  }

  const F_FLOAT a11 = -c*sb1*s1;
  const F_FLOAT a22 = sb2 * (2.0*costh13*s12 - c*(s1+s2));
  const F_FLOAT a33 = -c*sb3*s2;
  const F_FLOAT a12 = r12c1 * (costh12*c*s1 + costh23*s12);
  const F_FLOAT a13 = rb1*rb3*s12;
  const F_FLOAT a23 = r12c2 * (-costh23*c*s2 - costh12*s12);

  const F_FLOAT sx1  = a11*vb1x + a12*vb2x + a13*vb3x;
  const F_FLOAT sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
  const F_FLOAT sx12 = a13*vb1x + a23*vb2x + a33*vb3x;
  const F_FLOAT sy1  = a11*vb1y + a12*vb2y + a13*vb3y;
  const F_FLOAT sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
  const F_FLOAT sy12 = a13*vb1y + a23*vb2y + a33*vb3y;
  const F_FLOAT sz1  = a11*vb1z + a12*vb2z + a13*vb3z;
  const F_FLOAT sz2  = a12*vb1z + a22*vb2z + a23*vb3z;
  const F_FLOAT sz12 = a13*vb1z + a23*vb2z + a33*vb3z;

  // set up d(cos(phi))/d(r) and dphi/dr arrays

  F_FLOAT dcosphidr[4][3], dphidr[4][3];

  dcosphidr[0][0] = -sx1;
  dcosphidr[0][1] = -sy1;
  dcosphidr[0][2] = -sz1;
  dcosphidr[1][0] = sx2 + sx1;
  dcosphidr[1][1] = sy2 + sy1;
  dcosphidr[1][2] = sz2 + sz1;
  dcosphidr[2][0] = sx12 - sx2;
  dcosphidr[2][1] = sy12 - sy2;
  dcosphidr[2][2] = sz12 - sz2;
  dcosphidr[3][0] = -sx12;
  dcosphidr[3][1] = -sy12;
  dcosphidr[3][2] = -sz12;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      dphidr[i][j] = -dcosphidr[i][j] / sinphi;

  // energy

  F_FLOAT edihedral = 0.0;

  const F_FLOAT dphi1 = phi - d_phi1[type];
  const F_FLOAT dphi2 = 2.0*phi - d_phi2[type];
  const F_FLOAT dphi3 = 3.0*phi - d_phi3[type];

  if (eflag) edihedral = d_k1[type]*(1.0 - cos(dphi1)) +
                     d_k2[type]*(1.0 - cos(dphi2)) +
                     d_k3[type]*(1.0 - cos(dphi3));

  const F_FLOAT de_dihedral = d_k1[type]*sin(dphi1) + 2.0*d_k2[type]*sin(dphi2) +
      3.0*d_k3[type]*sin(dphi3);

  // torsion forces on all 4 atoms

  F_FLOAT dbonddr[3][4][3], fabcd[4][3];

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] = de_dihedral*dphidr[i][j];

  // set up d(bond)/d(r) array
  // dbonddr(i,j,k) = bond i, atom j, coordinate k

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 3; k++)
        dbonddr[i][j][k] = 0.0;

  // bond1

  dbonddr[0][0][0] = vb1x / r1;
  dbonddr[0][0][1] = vb1y / r1;
  dbonddr[0][0][2] = vb1z / r1;
  dbonddr[0][1][0] = -vb1x / r1;
  dbonddr[0][1][1] = -vb1y / r1;
  dbonddr[0][1][2] = -vb1z / r1;

  // bond2

  dbonddr[1][1][0] = vb2x / r2;
  dbonddr[1][1][1] = vb2y / r2;
  dbonddr[1][1][2] = vb2z / r2;
  dbonddr[1][2][0] = -vb2x / r2;
  dbonddr[1][2][1] = -vb2y / r2;
  dbonddr[1][2][2] = -vb2z / r2;

  // bond3

  dbonddr[2][2][0] = vb3x / r3;
  dbonddr[2][2][1] = vb3y / r3;
  dbonddr[2][2][2] = vb3z / r3;
  dbonddr[2][3][0] = -vb3x / r3;
  dbonddr[2][3][1] = -vb3y / r3;
  dbonddr[2][3][2] = -vb3z / r3;

  // set up d(theta)/d(r) array
  // dthetadr(i,j,k) = angle i, atom j, coordinate k

  F_FLOAT dthetadr[2][4][3];

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 3; k++)
        dthetadr[i][j][k] = 0.0;

  const F_FLOAT t1 = costh12 / r1mag2;
  const F_FLOAT t2 = costh23 / r2mag2;
  const F_FLOAT t3 = costh12 / r2mag2;
  const F_FLOAT t4 = costh23 / r3mag2;

  // angle12

  dthetadr[0][0][0] = sc1 * ((t1 * vb1x) - (vb2x * r12c1));
  dthetadr[0][0][1] = sc1 * ((t1 * vb1y) - (vb2y * r12c1));
  dthetadr[0][0][2] = sc1 * ((t1 * vb1z) - (vb2z * r12c1));

  dthetadr[0][1][0] = sc1 * ((-t1 * vb1x) + (vb2x * r12c1) +
                             (-t3 * vb2x) + (vb1x * r12c1));
  dthetadr[0][1][1] = sc1 * ((-t1 * vb1y) + (vb2y * r12c1) +
                 (-t3 * vb2y) + (vb1y * r12c1));
  dthetadr[0][1][2] = sc1 * ((-t1 * vb1z) + (vb2z * r12c1) +
                             (-t3 * vb2z) + (vb1z * r12c1));

  dthetadr[0][2][0] = sc1 * ((t3 * vb2x) - (vb1x * r12c1));
  dthetadr[0][2][1] = sc1 * ((t3 * vb2y) - (vb1y * r12c1));
  dthetadr[0][2][2] = sc1 * ((t3 * vb2z) - (vb1z * r12c1));

  // angle23

  dthetadr[1][1][0] = sc2 * ((t2 * vb2x) + (vb3x * r12c2));
  dthetadr[1][1][1] = sc2 * ((t2 * vb2y) + (vb3y * r12c2));
  dthetadr[1][1][2] = sc2 * ((t2 * vb2z) + (vb3z * r12c2));

  dthetadr[1][2][0] = sc2 * ((-t2 * vb2x) - (vb3x * r12c2) +
                          (t4 * vb3x) + (vb2x * r12c2));
  dthetadr[1][2][1] = sc2 * ((-t2 * vb2y) - (vb3y * r12c2) +
                          (t4 * vb3y) + (vb2y * r12c2));
  dthetadr[1][2][2] = sc2 * ((-t2 * vb2z) - (vb3z * r12c2) +
                          (t4 * vb3z) + (vb2z * r12c2));

  dthetadr[1][3][0] = -sc2 * ((t4 * vb3x) + (vb2x * r12c2));
  dthetadr[1][3][1] = -sc2 * ((t4 * vb3y) + (vb2y * r12c2));
  dthetadr[1][3][2] = -sc2 * ((t4 * vb3z) + (vb2z * r12c2));

  // mid-bond/torsion coupling
  // energy on bond2 (middle bond)

  F_FLOAT cos2phi = cos(2.0*phi);
  F_FLOAT cos3phi = cos(3.0*phi);

  F_FLOAT bt1 = d_mbt_f1[type] * cosphi;
  F_FLOAT bt2 = d_mbt_f2[type] * cos2phi;
  F_FLOAT bt3 = d_mbt_f3[type] * cos3phi;
  F_FLOAT sumbte = bt1 + bt2 + bt3;
  F_FLOAT db = r2 - d_mbt_r0[type];
  if (eflag) edihedral += db * sumbte;

  // force on bond2

  bt1 = -d_mbt_f1[type] * sinphi;
  bt2 = -2.0 * d_mbt_f2[type] * sin(2.0*phi);
  bt3 = -3.0 * d_mbt_f3[type] * sin(3.0*phi);
  F_FLOAT sumbtf = bt1 + bt2 + bt3;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] += db*sumbtf*dphidr[i][j] + sumbte*dbonddr[1][i][j];

  // end-bond/torsion coupling
  // energy on bond1 (first bond)

  bt1 = d_ebt_f1_1[type] * cosphi;
  bt2 = d_ebt_f2_1[type] * cos2phi;
  bt3 = d_ebt_f3_1[type] * cos3phi;
  sumbte = bt1 + bt2 + bt3;

  db = r1 - d_ebt_r0_1[type];
  if (eflag) edihedral += db * (bt1+bt2+bt3);

  // force on bond1

  bt1 = d_ebt_f1_1[type] * sinphi;
  bt2 = 2.0 * d_ebt_f2_1[type] * sin(2.0*phi);
  bt3 = 3.0 * d_ebt_f3_1[type] * sin(3.0*phi);
  sumbtf = bt1 + bt2 + bt3;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] -= db*sumbtf*dphidr[i][j] + sumbte*dbonddr[0][i][j];

  // end-bond/torsion coupling
  // energy on bond3 (last bond)

  bt1 = d_ebt_f1_2[type] * cosphi;
  bt2 = d_ebt_f2_2[type] * cos2phi;
  bt3 = d_ebt_f3_2[type] * cos3phi;
  sumbte = bt1 + bt2 + bt3;

  db = r3 - d_ebt_r0_2[type];
  if (eflag) edihedral += db * (bt1+bt2+bt3);

  // force on bond3

  bt1 = -d_ebt_f1_2[type] * sinphi;
  bt2 = -2.0 * d_ebt_f2_2[type] * sin(2.0*phi);
  bt3 = -3.0 * d_ebt_f3_2[type] * sin(3.0*phi);
  sumbtf = bt1 + bt2 + bt3;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] += db*sumbtf*dphidr[i][j] + sumbte*dbonddr[2][i][j];

  // angle/torsion coupling
  // energy on angle1

  F_FLOAT at1 = d_at_f1_1[type] * cosphi;
  F_FLOAT at2 = d_at_f2_1[type] * cos2phi;
  F_FLOAT at3 = d_at_f3_1[type] * cos3phi;
  sumbte = at1 + at2 + at3;

  F_FLOAT da = acos(costh12) - d_at_theta0_1[type];
  if (eflag) edihedral += da * (at1+at2+at3);

  // force on angle1

  bt1 = d_at_f1_1[type] * sinphi;
  bt2 = 2.0 * d_at_f2_1[type] * sin(2.0*phi);
  bt3 = 3.0 * d_at_f3_1[type] * sin(3.0*phi);
  sumbtf = bt1 + bt2 + bt3;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] -= da*sumbtf*dphidr[i][j] + sumbte*dthetadr[0][i][j];

  // energy on angle2

  at1 = d_at_f1_2[type] * cosphi;
  at2 = d_at_f2_2[type] * cos2phi;
  at3 = d_at_f3_2[type] * cos3phi;
  sumbte = at1 + at2 + at3;

  da = acos(costh23) - d_at_theta0_2[type];
  if (eflag) edihedral += da * (at1+at2+at3);

  // force on angle2

  bt1 = -d_at_f1_2[type] * sinphi;
  bt2 = -2.0 * d_at_f2_2[type] * sin(2.0*phi);
  bt3 = -3.0 * d_at_f3_2[type] * sin(3.0*phi);
  sumbtf = bt1 + bt2 + bt3;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] += da*sumbtf*dphidr[i][j] + sumbte*dthetadr[1][i][j];

  // angle/angle/torsion coupling

  const F_FLOAT da1 = acos(costh12) - d_aat_theta0_1[type];
  const F_FLOAT da2 = acos(costh23) - d_aat_theta0_2[type];

  if (eflag) edihedral += d_aat_k[type]*da1*da2*cosphi;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] -= d_aat_k[type] *
        (cosphi * (da2*dthetadr[0][i][j] - da1*dthetadr[1][i][j]) +
         sinphi * da1*da2*dphidr[i][j]);

  // bond1/bond3 coupling

  if (fabs(d_bb13t_k[type]) > SMALL) {

    const F_FLOAT r1_0 = d_bb13t_r10[type];
    const F_FLOAT r3_0 = d_bb13t_r30[type];
    const F_FLOAT dr1 = r1 - r1_0;
    const F_FLOAT dr2 = r3 - r3_0;
    const F_FLOAT tk1 = -d_bb13t_k[type] * dr1 / r3;
    const F_FLOAT tk2 = -d_bb13t_k[type] * dr2 / r1;

    if (eflag) edihedral += d_bb13t_k[type]*dr1*dr2;

    fabcd[0][0] += tk2 * vb1x;
    fabcd[0][1] += tk2 * vb1y;
    fabcd[0][2] += tk2 * vb1z;

    fabcd[1][0] -= tk2 * vb1x;
    fabcd[1][1] -= tk2 * vb1y;
    fabcd[1][2] -= tk2 * vb1z;

    fabcd[2][0] -= tk1 * vb3x;
    fabcd[2][1] -= tk1 * vb3y;
    fabcd[2][2] -= tk1 * vb3z;

    fabcd[3][0] += tk1 * vb3x;
    fabcd[3][1] += tk1 * vb3y;
    fabcd[3][2] += tk1 * vb3z;
  }

  F_FLOAT f1[3],f2[3],f3[3],f4[3];

  for (int i = 0; i < 3; i++) {
    f1[i] = fabcd[0][i];
    f2[i] = fabcd[1][i];
    f3[i] = fabcd[2][i];
    f4[i] = fabcd[3][i];
  }

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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void DihedralClass2Kokkos<DeviceType>::operator()(TagDihedralClass2Compute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralClass2Compute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralClass2Kokkos<DeviceType>::allocate()
{
  DihedralClass2::allocate();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralClass2Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  DihedralClass2::coeff(narg, arg);

  int n = atom->ndihedraltypes;
  k_k1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::k1",n+1);
  k_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::k2",n+1);
  k_k3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::k3",n+1);
  k_phi1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::phi1",n+1);
  k_phi2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::phi2",n+1);
  k_phi3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::phi3",n+1);
  k_mbt_f1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_f1",n+1);
  k_mbt_f2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_f2",n+1);
  k_mbt_f3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_f3",n+1);
  k_mbt_r0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_r0",n+1);
  k_ebt_f1_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f1_1",n+1);
  k_ebt_f2_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f2_1",n+1);
  k_ebt_f3_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f3_1",n+1);
  k_ebt_r0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_r0_1",n+1);
  k_ebt_f1_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f1_2",n+1);
  k_ebt_f2_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f2_2",n+1);
  k_ebt_f3_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f3_2",n+1);
  k_ebt_r0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_r0_2",n+1);
  k_at_f1_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f1_1",n+1);
  k_at_f2_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f2_1",n+1);
  k_at_f3_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f3_1",n+1);
  k_at_f1_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f1_2",n+1);
  k_at_f2_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f2_2",n+1);
  k_at_f3_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f3_2",n+1);
  k_at_theta0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_theta0_1",n+1);
  k_at_theta0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_theta0_2",n+1);
  k_aat_k = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::aat_k",n+1);
  k_aat_theta0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::aat_theta0_1",n+1);
  k_aat_theta0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::aat_theta0_2",n+1);
  k_bb13t_k = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::bb13t_k",n+1);
  k_bb13t_r10 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::bb13t_r10",n+1);
  k_bb13t_r30 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::bb13t_r30",n+1);
  k_setflag_d = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_d",n+1);
  k_setflag_mbt = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_mbt",n+1);
  k_setflag_ebt = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_ebt",n+1);
  k_setflag_at = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_at",n+1);
  k_setflag_aat = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_aat",n+1);
  k_setflag_bb13t = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_bb13t",n+1);

  d_k1 = k_k1.template view<DeviceType>();
  d_k2 = k_k2.template view<DeviceType>();
  d_k3 = k_k3.template view<DeviceType>();
  d_phi1 = k_phi1.template view<DeviceType>();
  d_phi2 = k_phi2.template view<DeviceType>();
  d_phi3 = k_phi3.template view<DeviceType>();
  d_mbt_f1 = k_mbt_f1.template view<DeviceType>();
  d_mbt_f2 = k_mbt_f2.template view<DeviceType>();
  d_mbt_f3 = k_mbt_f3.template view<DeviceType>();
  d_mbt_r0 = k_mbt_r0.template view<DeviceType>();
  d_ebt_f1_1 = k_ebt_f1_1.template view<DeviceType>();
  d_ebt_f2_1 = k_ebt_f2_1.template view<DeviceType>();
  d_ebt_f3_1 = k_ebt_f3_1.template view<DeviceType>();
  d_ebt_r0_1 = k_ebt_r0_1.template view<DeviceType>();
  d_ebt_f1_2 = k_ebt_f1_2.template view<DeviceType>();
  d_ebt_f2_2 = k_ebt_f2_2.template view<DeviceType>();
  d_ebt_f3_2 = k_ebt_f3_2.template view<DeviceType>();
  d_ebt_r0_2 = k_ebt_r0_2.template view<DeviceType>();
  d_at_f1_1 = k_at_f1_1.template view<DeviceType>();
  d_at_f2_1 = k_at_f2_1.template view<DeviceType>();
  d_at_f3_1 = k_at_f3_1.template view<DeviceType>();
  d_at_f1_2 = k_at_f1_2.template view<DeviceType>();
  d_at_f2_2 = k_at_f2_2.template view<DeviceType>();
  d_at_f3_2 = k_at_f3_2.template view<DeviceType>();
  d_at_theta0_1 = k_at_theta0_1.template view<DeviceType>();
  d_at_theta0_2 = k_at_theta0_2.template view<DeviceType>();
  d_aat_k = k_aat_k.template view<DeviceType>();
  d_aat_theta0_1 = k_aat_theta0_1.template view<DeviceType>();
  d_aat_theta0_2 = k_aat_theta0_2.template view<DeviceType>();
  d_bb13t_k = k_bb13t_k.template view<DeviceType>();
  d_bb13t_r10 = k_bb13t_r10.template view<DeviceType>();
  d_bb13t_r30 = k_bb13t_r30.template view<DeviceType>();
  d_setflag_d = k_setflag_d.template view<DeviceType>();
  d_setflag_mbt = k_setflag_mbt.template view<DeviceType>();
  d_setflag_ebt = k_setflag_ebt.template view<DeviceType>();
  d_setflag_at = k_setflag_at.template view<DeviceType>();
  d_setflag_aat = k_setflag_aat.template view<DeviceType>();
  d_setflag_bb13t = k_setflag_bb13t.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k1.h_view[i] = k1[i];
    k_k2.h_view[i] = k2[i];
    k_k3.h_view[i] = k3[i];
    k_phi1.h_view[i] = phi1[i];
    k_phi2.h_view[i] = phi2[i];
    k_phi3.h_view[i] = phi3[i];
    k_mbt_f1.h_view[i] = mbt_f1[i];
    k_mbt_f2.h_view[i] = mbt_f2[i];
    k_mbt_f3.h_view[i] = mbt_f3[i];
    k_mbt_r0.h_view[i] = mbt_r0[i];
    k_ebt_f1_1.h_view[i] = ebt_f1_1[i];
    k_ebt_f2_1.h_view[i] = ebt_f2_1[i];
    k_ebt_f3_1.h_view[i] = ebt_f3_1[i];
    k_ebt_r0_1.h_view[i] = ebt_r0_1[i];
    k_ebt_f1_2.h_view[i] = ebt_f1_2[i];
    k_ebt_f2_2.h_view[i] = ebt_f2_2[i];
    k_ebt_f3_2.h_view[i] = ebt_f3_2[i];
    k_ebt_r0_2.h_view[i] = ebt_r0_2[i];
    k_at_f1_1.h_view[i] = at_f1_1[i];
    k_at_f2_1.h_view[i] = at_f2_1[i];
    k_at_f3_1.h_view[i] = at_f3_1[i];
    k_at_f1_2.h_view[i] = at_f1_2[i];
    k_at_f2_2.h_view[i] = at_f2_2[i];
    k_at_f3_2.h_view[i] = at_f3_2[i];
    k_at_theta0_1.h_view[i] = at_theta0_1[i];
    k_at_theta0_2.h_view[i] = at_theta0_2[i];
    k_aat_k.h_view[i] = aat_k[i];
    k_aat_theta0_1.h_view[i] = aat_theta0_1[i];
    k_aat_theta0_2.h_view[i] = aat_theta0_2[i];
    k_bb13t_k.h_view[i] = bb13t_k[i];
    k_bb13t_r10.h_view[i] = bb13t_r10[i];
    k_bb13t_r30.h_view[i] = bb13t_r30[i];
    k_setflag_d.h_view[i] = setflag_d[i];
    k_setflag_mbt.h_view[i] = setflag_mbt[i];
    k_setflag_ebt.h_view[i] = setflag_ebt[i];
    k_setflag_at.h_view[i] = setflag_at[i];
    k_setflag_aat.h_view[i] = setflag_aat[i];
    k_setflag_bb13t.h_view[i] = setflag_bb13t[i];
  }

  k_k1.template modify<LMPHostType>();
  k_k2.template modify<LMPHostType>();
  k_k3.template modify<LMPHostType>();
  k_phi1.template modify<LMPHostType>();
  k_phi2.template modify<LMPHostType>();
  k_phi3.template modify<LMPHostType>();
  k_mbt_f1.template modify<LMPHostType>();
  k_mbt_f2.template modify<LMPHostType>();
  k_mbt_f3.template modify<LMPHostType>();
  k_mbt_r0.template modify<LMPHostType>();
  k_ebt_f1_1.template modify<LMPHostType>();
  k_ebt_f2_1.template modify<LMPHostType>();
  k_ebt_f3_1.template modify<LMPHostType>();
  k_ebt_r0_1.template modify<LMPHostType>();
  k_ebt_f1_2.template modify<LMPHostType>();
  k_ebt_f2_2.template modify<LMPHostType>();
  k_ebt_f3_2.template modify<LMPHostType>();
  k_ebt_r0_2.template modify<LMPHostType>();
  k_at_f1_1.template modify<LMPHostType>();
  k_at_f2_1.template modify<LMPHostType>();
  k_at_f3_1.template modify<LMPHostType>();
  k_at_f1_2.template modify<LMPHostType>();
  k_at_f2_2.template modify<LMPHostType>();
  k_at_f3_2.template modify<LMPHostType>();
  k_at_theta0_1.template modify<LMPHostType>();
  k_at_theta0_2.template modify<LMPHostType>();
  k_aat_k.template modify<LMPHostType>();
  k_aat_theta0_1.template modify<LMPHostType>();
  k_aat_theta0_2.template modify<LMPHostType>();
  k_bb13t_k.template modify<LMPHostType>();
  k_bb13t_r10.template modify<LMPHostType>();
  k_bb13t_r30.template modify<LMPHostType>();
  k_setflag_d.template modify<LMPHostType>();
  k_setflag_mbt.template modify<LMPHostType>();
  k_setflag_ebt.template modify<LMPHostType>();
  k_setflag_at.template modify<LMPHostType>();
  k_setflag_aat.template modify<LMPHostType>();
  k_setflag_bb13t.template modify<LMPHostType>();
}


/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralClass2Kokkos<DeviceType>::read_restart(FILE *fp)
{
  DihedralClass2::read_restart(fp);

  int n = atom->ndihedraltypes;
  k_k1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::k1",n+1);
  k_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::k2",n+1);
  k_k3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::k3",n+1);
  k_phi1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::phi1",n+1);
  k_phi2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::phi2",n+1);
  k_phi3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::phi3",n+1);
  k_mbt_f1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_f1",n+1);
  k_mbt_f2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_f2",n+1);
  k_mbt_f3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_f3",n+1);
  k_mbt_r0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::mbt_r0",n+1);
  k_ebt_f1_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f1_1",n+1);
  k_ebt_f2_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f2_1",n+1);
  k_ebt_f3_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f3_1",n+1);
  k_ebt_r0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_r0_1",n+1);
  k_ebt_f1_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f1_2",n+1);
  k_ebt_f2_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f2_2",n+1);
  k_ebt_f3_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_f3_2",n+1);
  k_ebt_r0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::ebt_r0_2",n+1);
  k_at_f1_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f1_1",n+1);
  k_at_f2_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f2_1",n+1);
  k_at_f3_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f3_1",n+1);
  k_at_f1_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f1_2",n+1);
  k_at_f2_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f2_2",n+1);
  k_at_f3_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_f3_2",n+1);
  k_at_theta0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_theta0_1",n+1);
  k_at_theta0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::at_theta0_2",n+1);
  k_aat_k = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::aat_k",n+1);
  k_aat_theta0_1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::aat_theta0_1",n+1);
  k_aat_theta0_2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::aat_theta0_2",n+1);
  k_bb13t_k = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::bb13t_k",n+1);
  k_bb13t_r10 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::bb13t_r10",n+1);
  k_bb13t_r30 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("DihedralClass2::bb13t_r30",n+1);
  k_setflag_d = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_d",n+1);
  k_setflag_mbt = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_mbt",n+1);
  k_setflag_ebt = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_ebt",n+1);
  k_setflag_at = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_at",n+1);
  k_setflag_aat = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_aat",n+1);
  k_setflag_bb13t = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_bb13t",n+1);

  d_k1 = k_k1.template view<DeviceType>();
  d_k2 = k_k2.template view<DeviceType>();
  d_k3 = k_k3.template view<DeviceType>();
  d_phi1 = k_phi1.template view<DeviceType>();
  d_phi2 = k_phi2.template view<DeviceType>();
  d_phi3 = k_phi3.template view<DeviceType>();
  d_mbt_f1 = k_mbt_f1.template view<DeviceType>();
  d_mbt_f2 = k_mbt_f2.template view<DeviceType>();
  d_mbt_f3 = k_mbt_f3.template view<DeviceType>();
  d_mbt_r0 = k_mbt_r0.template view<DeviceType>();
  d_ebt_f1_1 = k_ebt_f1_1.template view<DeviceType>();
  d_ebt_f2_1 = k_ebt_f2_1.template view<DeviceType>();
  d_ebt_f3_1 = k_ebt_f3_1.template view<DeviceType>();
  d_ebt_r0_1 = k_ebt_r0_1.template view<DeviceType>();
  d_ebt_f1_2 = k_ebt_f1_2.template view<DeviceType>();
  d_ebt_f2_2 = k_ebt_f2_2.template view<DeviceType>();
  d_ebt_f3_2 = k_ebt_f3_2.template view<DeviceType>();
  d_ebt_r0_2 = k_ebt_r0_2.template view<DeviceType>();
  d_at_f1_1 = k_at_f1_1.template view<DeviceType>();
  d_at_f2_1 = k_at_f2_1.template view<DeviceType>();
  d_at_f3_1 = k_at_f3_1.template view<DeviceType>();
  d_at_f1_2 = k_at_f1_2.template view<DeviceType>();
  d_at_f2_2 = k_at_f2_2.template view<DeviceType>();
  d_at_f3_2 = k_at_f3_2.template view<DeviceType>();
  d_at_theta0_1 = k_at_theta0_1.template view<DeviceType>();
  d_at_theta0_2 = k_at_theta0_2.template view<DeviceType>();
  d_aat_k = k_aat_k.template view<DeviceType>();
  d_aat_theta0_1 = k_aat_theta0_1.template view<DeviceType>();
  d_aat_theta0_2 = k_aat_theta0_2.template view<DeviceType>();
  d_bb13t_k = k_bb13t_k.template view<DeviceType>();
  d_bb13t_r10 = k_bb13t_r10.template view<DeviceType>();
  d_bb13t_r30 = k_bb13t_r30.template view<DeviceType>();
  d_setflag_d = k_setflag_d.template view<DeviceType>();
  d_setflag_mbt = k_setflag_mbt.template view<DeviceType>();
  d_setflag_ebt = k_setflag_ebt.template view<DeviceType>();
  d_setflag_at = k_setflag_at.template view<DeviceType>();
  d_setflag_aat = k_setflag_aat.template view<DeviceType>();
  d_setflag_bb13t = k_setflag_bb13t.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k1.h_view[i] = k1[i];
    k_k2.h_view[i] = k2[i];
    k_k3.h_view[i] = k3[i];
    k_phi1.h_view[i] = phi1[i];
    k_phi2.h_view[i] = phi2[i];
    k_phi3.h_view[i] = phi3[i];
    k_mbt_f1.h_view[i] = mbt_f1[i];
    k_mbt_f2.h_view[i] = mbt_f2[i];
    k_mbt_f3.h_view[i] = mbt_f3[i];
    k_mbt_r0.h_view[i] = mbt_r0[i];
    k_ebt_f1_1.h_view[i] = ebt_f1_1[i];
    k_ebt_f2_1.h_view[i] = ebt_f2_1[i];
    k_ebt_f3_1.h_view[i] = ebt_f3_1[i];
    k_ebt_r0_1.h_view[i] = ebt_r0_1[i];
    k_ebt_f1_2.h_view[i] = ebt_f1_2[i];
    k_ebt_f2_2.h_view[i] = ebt_f2_2[i];
    k_ebt_f3_2.h_view[i] = ebt_f3_2[i];
    k_ebt_r0_2.h_view[i] = ebt_r0_2[i];
    k_at_f1_1.h_view[i] = at_f1_1[i];
    k_at_f2_1.h_view[i] = at_f2_1[i];
    k_at_f3_1.h_view[i] = at_f3_1[i];
    k_at_f1_2.h_view[i] = at_f1_2[i];
    k_at_f2_2.h_view[i] = at_f2_2[i];
    k_at_f3_2.h_view[i] = at_f3_2[i];
    k_at_theta0_1.h_view[i] = at_theta0_1[i];
    k_at_theta0_2.h_view[i] = at_theta0_2[i];
    k_aat_k.h_view[i] = aat_k[i];
    k_aat_theta0_1.h_view[i] = aat_theta0_1[i];
    k_aat_theta0_2.h_view[i] = aat_theta0_2[i];
    k_bb13t_k.h_view[i] = bb13t_k[i];
    k_bb13t_r10.h_view[i] = bb13t_r10[i];
    k_bb13t_r30.h_view[i] = bb13t_r30[i];
    k_setflag_d.h_view[i] = setflag_d[i];
    k_setflag_mbt.h_view[i] = setflag_mbt[i];
    k_setflag_ebt.h_view[i] = setflag_ebt[i];
    k_setflag_at.h_view[i] = setflag_at[i];
    k_setflag_aat.h_view[i] = setflag_aat[i];
    k_setflag_bb13t.h_view[i] = setflag_bb13t[i];
  }

  k_k1.template modify<LMPHostType>();
  k_k2.template modify<LMPHostType>();
  k_k3.template modify<LMPHostType>();
  k_phi1.template modify<LMPHostType>();
  k_phi2.template modify<LMPHostType>();
  k_phi3.template modify<LMPHostType>();
  k_mbt_f1.template modify<LMPHostType>();
  k_mbt_f2.template modify<LMPHostType>();
  k_mbt_f3.template modify<LMPHostType>();
  k_mbt_r0.template modify<LMPHostType>();
  k_ebt_f1_1.template modify<LMPHostType>();
  k_ebt_f2_1.template modify<LMPHostType>();
  k_ebt_f3_1.template modify<LMPHostType>();
  k_ebt_r0_1.template modify<LMPHostType>();
  k_ebt_f1_2.template modify<LMPHostType>();
  k_ebt_f2_2.template modify<LMPHostType>();
  k_ebt_f3_2.template modify<LMPHostType>();
  k_ebt_r0_2.template modify<LMPHostType>();
  k_at_f1_1.template modify<LMPHostType>();
  k_at_f2_1.template modify<LMPHostType>();
  k_at_f3_1.template modify<LMPHostType>();
  k_at_f1_2.template modify<LMPHostType>();
  k_at_f2_2.template modify<LMPHostType>();
  k_at_f3_2.template modify<LMPHostType>();
  k_at_theta0_1.template modify<LMPHostType>();
  k_at_theta0_2.template modify<LMPHostType>();
  k_aat_k.template modify<LMPHostType>();
  k_aat_theta0_1.template modify<LMPHostType>();
  k_aat_theta0_2.template modify<LMPHostType>();
  k_bb13t_k.template modify<LMPHostType>();
  k_bb13t_r10.template modify<LMPHostType>();
  k_bb13t_r30.template modify<LMPHostType>();
  k_setflag_d.template modify<LMPHostType>();
  k_setflag_mbt.template modify<LMPHostType>();
  k_setflag_ebt.template modify<LMPHostType>();
  k_setflag_at.template modify<LMPHostType>();
  k_setflag_aat.template modify<LMPHostType>();
  k_setflag_bb13t.template modify<LMPHostType>();
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
void DihedralClass2Kokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                        F_FLOAT &edihedral, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                        const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                        const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                        const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const
{
  E_FLOAT edihedralquarter;
  F_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.view<DeviceType>();

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
template class DihedralClass2Kokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class DihedralClass2Kokkos<LMPHostType>;
#endif
}
