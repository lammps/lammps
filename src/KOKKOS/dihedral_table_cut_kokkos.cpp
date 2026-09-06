// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "dihedral_table_cut_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int LINEAR_STYLE = 0;
static constexpr double TOLERANCE = 0.05;
static constexpr double SMALL = 0.0000001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralTableCutKokkos<DeviceType>::DihedralTableCutKokkos(LAMMPS *lmp) : DihedralTableCut(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = DAT::tdual_int_scalar("Dihedral:warning_flag");
  d_warning_flag = k_warning_flag.template view<DeviceType>();
  h_warning_flag = k_warning_flag.view_host();


  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralTableCutKokkos<DeviceType>::~DihedralTableCutKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralTableCutKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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


  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_dihedrallist.template sync<DeviceType>();
  dihedrallist = neighborKK->k_dihedrallist.view<DeviceType>();
  int ndihedrallist = neighborKK->ndihedrallist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  h_warning_flag() = 0;
  k_warning_flag.modify_host();
  k_warning_flag.template sync<DeviceType>();

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralTableCutCompute<1,1> >(0,ndihedrallist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralTableCutCompute<0,1> >(0,ndihedrallist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralTableCutCompute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralTableCutCompute<0,0> >(0,ndihedrallist),*this);
    }
  }

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.sync_host();
  if (h_warning_flag())
    error->warning(FLERR,"Dihedral problem");

  if (eflag_global) energy += static_cast<double>(ev.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  copymode = 0;
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralTableCutKokkos<DeviceType>::operator()(TagDihedralTableCutCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {
  // The f array is atomic
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = dihedrallist(n,0);
  const int i2 = dihedrallist(n,1);
  const int i3 = dihedrallist(n,2);
  const int i4 = dihedrallist(n,3);
  const int type = dihedrallist(n,4);

  // 1st bond

  const KK_FLOAT vb1x = x(i1,0) - x(i2,0);
  const KK_FLOAT vb1y = x(i1,1) - x(i2,1);
  const KK_FLOAT vb1z = x(i1,2) - x(i2,2);

  // 2nd bond

  const KK_FLOAT vb2x = x(i3,0) - x(i2,0);
  const KK_FLOAT vb2y = x(i3,1) - x(i2,1);
  const KK_FLOAT vb2z = x(i3,2) - x(i2,2);

  const KK_FLOAT vb2xm = -vb2x;
  const KK_FLOAT vb2ym = -vb2y;
  const KK_FLOAT vb2zm = -vb2z;

  // 3rd bond

  const KK_FLOAT vb3x = x(i4,0) - x(i3,0);
  const KK_FLOAT vb3y = x(i4,1) - x(i3,1);
  const KK_FLOAT vb3z = x(i4,2) - x(i3,2);

  // distances

  const KK_FLOAT r1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
  const KK_FLOAT r1 = Kokkos::sqrt(r1mag2);
  const KK_FLOAT r2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
  const KK_FLOAT r2 = Kokkos::sqrt(r2mag2);
  const KK_FLOAT r3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
  const KK_FLOAT r3 = Kokkos::sqrt(r3mag2);

  const KK_FLOAT sb1 = static_cast<KK_FLOAT>(1.0)/r1mag2;
  const KK_FLOAT rb1 = static_cast<KK_FLOAT>(1.0)/r1;
  const KK_FLOAT sb2 = static_cast<KK_FLOAT>(1.0)/r2mag2;
  const KK_FLOAT rb2 = static_cast<KK_FLOAT>(1.0)/r2;
  const KK_FLOAT sb3 = static_cast<KK_FLOAT>(1.0)/r3mag2;
  const KK_FLOAT rb3 = static_cast<KK_FLOAT>(1.0)/r3;

  const KK_FLOAT c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

  // angles

  const KK_FLOAT r12c1 = rb1*rb2;
  const KK_FLOAT r12c2 = rb2*rb3;
  const KK_FLOAT costh12 = (vb1x*vb2x + vb1y*vb2y + vb1z*vb2z) * r12c1;
  const KK_FLOAT costh13 = c0;
  const KK_FLOAT costh23 = (vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z) * r12c2;

  // cos and sin of 2 angles and final c

  KK_FLOAT sin2 = Kokkos::fmax(static_cast<KK_FLOAT>(1.0) - costh12*costh12,static_cast<KK_FLOAT>(0.0));
  KK_FLOAT sc1 = Kokkos::sqrt(sin2);
  if (sc1 < static_cast<KK_FLOAT>(SMALL)) sc1 = static_cast<KK_FLOAT>(SMALL);
  sc1 = static_cast<KK_FLOAT>(1.0)/sc1;

  sin2 = Kokkos::fmax(static_cast<KK_FLOAT>(1.0) - costh23*costh23,static_cast<KK_FLOAT>(0.0));
  KK_FLOAT sc2 = Kokkos::sqrt(sin2);
  if (sc2 < static_cast<KK_FLOAT>(SMALL)) sc2 = static_cast<KK_FLOAT>(SMALL);
  sc2 = static_cast<KK_FLOAT>(1.0)/sc2;

  const KK_FLOAT s1 = sc1 * sc1;
  const KK_FLOAT s2 = sc2 * sc2;
  const KK_FLOAT s12 = sc1 * sc2;
  KK_FLOAT c = (c0 + costh12*costh23) * s12;

  // error check

  if ((c > static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(TOLERANCE) ||
       c < static_cast<KK_FLOAT>(-1.0) - static_cast<KK_FLOAT>(TOLERANCE)) && !d_warning_flag())
    d_warning_flag() = 1;

  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);
  KK_FLOAT phil = Kokkos::acos(c);

  KK_FLOAT sinphi = Kokkos::sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  sinphi = Kokkos::fmax(sinphi,static_cast<KK_FLOAT>(SMALL));

  // n123 = vb1 x vb2

  const KK_FLOAT n123x = vb1y*vb2z - vb1z*vb2y;
  const KK_FLOAT n123y = vb1z*vb2x - vb1x*vb2z;
  const KK_FLOAT n123z = vb1x*vb2y - vb1y*vb2x;
  const KK_FLOAT n123_dot_vb3 = n123x*vb3x + n123y*vb3y + n123z*vb3z;
  if (n123_dot_vb3 > static_cast<KK_FLOAT>(0.0)) {
    phil = -phil;
    sinphi = -sinphi;
  }

  const KK_FLOAT a11 = -c*sb1*s1;
  const KK_FLOAT a22 = sb2 * (static_cast<KK_FLOAT>(2.0)*costh13*s12 - c*(s1+s2));
  const KK_FLOAT a33 = -c*sb3*s2;
  const KK_FLOAT a12 = r12c1 * (costh12*c*s1 + costh23*s12);
  const KK_FLOAT a13 = rb1*rb3*s12;
  const KK_FLOAT a23 = r12c2 * (-costh23*c*s2 - costh12*s12);

  const KK_FLOAT sx1  = a11*vb1x + a12*vb2x + a13*vb3x;
  const KK_FLOAT sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
  const KK_FLOAT sx12 = a13*vb1x + a23*vb2x + a33*vb3x;
  const KK_FLOAT sy1  = a11*vb1y + a12*vb2y + a13*vb3y;
  const KK_FLOAT sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
  const KK_FLOAT sy12 = a13*vb1y + a23*vb2y + a33*vb3y;
  const KK_FLOAT sz1  = a11*vb1z + a12*vb2z + a13*vb3z;
  const KK_FLOAT sz2  = a12*vb1z + a22*vb2z + a23*vb3z;
  const KK_FLOAT sz12 = a13*vb1z + a23*vb2z + a33*vb3z;

  // set up d(cos(phi))/d(r) and dphi/dr arrays

  KK_FLOAT dcosphidr[4][3],dphidr[4][3],dthetadr[2][4][3],fabcd[4][3];

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
    for (int j = 0; j < 3; j++) {
      dphidr[i][j] = -dcosphidr[i][j] / sinphi;
      fabcd[i][j] = 0;
    }

  // set up d(theta)/d(r) array
  // dthetadr(i,j,k) = angle i, atom j, coordinate k

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 3; k++)
        dthetadr[i][j][k] = static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT t1 = costh12 / r1mag2;
  const KK_FLOAT t2 = costh23 / r2mag2;
  const KK_FLOAT t3 = costh12 / r2mag2;
  const KK_FLOAT t4 = costh23 / r3mag2;

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

  // angle/angle/torsion cutoff

  const KK_FLOAT aat_k = d_aat_k[type];
  const KK_FLOAT aat_theta0_1 = d_aat_theta0_1[type];
  const KK_FLOAT da1 = Kokkos::acos(costh12) - aat_theta0_1;
  const KK_FLOAT da2 = Kokkos::acos(costh23) - aat_theta0_1;
  const KK_FLOAT dtheta = d_aat_theta0_2[type] - aat_theta0_1;

  KK_FLOAT fphi = static_cast<KK_FLOAT>(0.0);
  KK_FLOAT fpphi = static_cast<KK_FLOAT>(0.0);
  if (phil < static_cast<KK_FLOAT>(0.0)) phil += static_cast<KK_FLOAT>(MY_2PI);
  uf_lookup_kk(type,phil,fphi,fpphi);

  KK_FLOAT gt = aat_k;
  KK_FLOAT gtt = aat_k;
  KK_FLOAT gpt = static_cast<KK_FLOAT>(0.0);
  KK_FLOAT gptt = static_cast<KK_FLOAT>(0.0);

  if (Kokkos::acos(costh12) > aat_theta0_1) {
    gt *= static_cast<KK_FLOAT>(1.0) - da1*da1/dtheta/dtheta;
    gpt = -aat_k*static_cast<KK_FLOAT>(2.0)*da1/dtheta/dtheta;
  }

  if (Kokkos::acos(costh23) > aat_theta0_1) {
    gtt *= static_cast<KK_FLOAT>(1.0) - da2*da2/dtheta/dtheta;
    gptt = -aat_k*static_cast<KK_FLOAT>(2.0)*da2/dtheta/dtheta;
  }

  KK_FLOAT edihedral = 0;
  if (eflag) edihedral = gt*gtt*fphi;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      fabcd[i][j] -= gt*gtt*fpphi*dphidr[i][j]
        - gt*gptt*fphi*dthetadr[1][i][j] + gpt*gtt*fphi*dthetadr[0][i][j];

  // apply force to each of 4 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += static_cast<KK_ACC_FLOAT>(fabcd[0][0]);
    a_f(i1,1) += static_cast<KK_ACC_FLOAT>(fabcd[0][1]);
    a_f(i1,2) += static_cast<KK_ACC_FLOAT>(fabcd[0][2]);
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) += static_cast<KK_ACC_FLOAT>(fabcd[1][0]);
    a_f(i2,1) += static_cast<KK_ACC_FLOAT>(fabcd[1][1]);
    a_f(i2,2) += static_cast<KK_ACC_FLOAT>(fabcd[1][2]);
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += static_cast<KK_ACC_FLOAT>(fabcd[2][0]);
    a_f(i3,1) += static_cast<KK_ACC_FLOAT>(fabcd[2][1]);
    a_f(i3,2) += static_cast<KK_ACC_FLOAT>(fabcd[2][2]);
  }

  if (NEWTON_BOND || i4 < nlocal) {
    a_f(i4,0) += static_cast<KK_ACC_FLOAT>(fabcd[3][0]);
    a_f(i4,1) += static_cast<KK_ACC_FLOAT>(fabcd[3][1]);
    a_f(i4,2) += static_cast<KK_ACC_FLOAT>(fabcd[3][2]);
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,edihedral,fabcd[0],fabcd[2],fabcd[3],
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralTableCutKokkos<DeviceType>::operator()(TagDihedralTableCutCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralTableCutCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ----------------------------------------------------------------------
   copy the tabulated data built on the host by compute_table() into
   device views.  this runs once per setup, after every coeff() call, so
   restarts and repeated dihedral_coeff commands are covered
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralTableCutKokkos<DeviceType>::setup_tables()
{
  const int n = atom->ndihedraltypes;

  k_tabindex = DAT::tdual_int_1d("DihedralTableCut::tabindex",n+1);
  k_f_unspecified = DAT::tdual_int_1d("DihedralTableCut::f_unspecified",ntables);
  k_delta = DAT::tdual_kkfloat_1d("DihedralTableCut::delta",ntables);
  k_invdelta = DAT::tdual_kkfloat_1d("DihedralTableCut::invdelta",ntables);
  k_deltasq6 = DAT::tdual_kkfloat_1d("DihedralTableCut::deltasq6",ntables);

  k_e = DAT::tdual_kkfloat_2d("DihedralTableCut::e",ntables,tablength);
  k_de = DAT::tdual_kkfloat_2d("DihedralTableCut::de",ntables,tablength);
  k_f_tab = DAT::tdual_kkfloat_2d("DihedralTableCut::f",ntables,tablength);
  k_df = DAT::tdual_kkfloat_2d("DihedralTableCut::df",ntables,tablength);
  k_e2 = DAT::tdual_kkfloat_2d("DihedralTableCut::e2",ntables,tablength);
  k_f2 = DAT::tdual_kkfloat_2d("DihedralTableCut::f2",ntables,tablength);

  k_aat_k = DAT::tdual_kkfloat_1d("DihedralTableCut::aat_k",n+1);
  k_aat_theta0_1 = DAT::tdual_kkfloat_1d("DihedralTableCut::aat_theta0_1",n+1);
  k_aat_theta0_2 = DAT::tdual_kkfloat_1d("DihedralTableCut::aat_theta0_2",n+1);

  for (int i = 1; i <= n; i++) {
    k_tabindex.view_host()(i) = tabindex[i];
    k_aat_k.view_host()(i) = static_cast<KK_FLOAT>(aat_k[i]);
    k_aat_theta0_1.view_host()(i) = static_cast<KK_FLOAT>(aat_theta0_1[i]);
    k_aat_theta0_2.view_host()(i) = static_cast<KK_FLOAT>(aat_theta0_2[i]);
  }

  for (int m = 0; m < ntables; m++) {
    const Table *tb = &tables[m];
    k_f_unspecified.view_host()(m) = tb->f_unspecified;
    k_delta.view_host()(m) = static_cast<KK_FLOAT>(tb->delta);
    k_invdelta.view_host()(m) = static_cast<KK_FLOAT>(tb->invdelta);
    k_deltasq6.view_host()(m) = static_cast<KK_FLOAT>(tb->deltasq6);
    for (int i = 0; i < tablength; i++) {
      k_e.view_host()(m,i) = static_cast<KK_FLOAT>(tb->e[i]);
      k_de.view_host()(m,i) = static_cast<KK_FLOAT>(tb->de[i]);
      k_f_tab.view_host()(m,i) = static_cast<KK_FLOAT>(tb->f[i]);
      k_df.view_host()(m,i) = static_cast<KK_FLOAT>(tb->df[i]);
      k_e2.view_host()(m,i) = static_cast<KK_FLOAT>(tb->e2[i]);
      k_f2.view_host()(m,i) = static_cast<KK_FLOAT>(tb->f2[i]);
    }
  }

  k_tabindex.modify_host(); k_tabindex.template sync<DeviceType>();
  k_aat_k.modify_host(); k_aat_k.template sync<DeviceType>();
  k_aat_theta0_1.modify_host(); k_aat_theta0_1.template sync<DeviceType>();
  k_aat_theta0_2.modify_host(); k_aat_theta0_2.template sync<DeviceType>();
  k_f_unspecified.modify_host(); k_f_unspecified.template sync<DeviceType>();
  k_delta.modify_host(); k_delta.template sync<DeviceType>();
  k_invdelta.modify_host(); k_invdelta.template sync<DeviceType>();
  k_deltasq6.modify_host(); k_deltasq6.template sync<DeviceType>();
  k_e.modify_host(); k_e.template sync<DeviceType>();
  k_de.modify_host(); k_de.template sync<DeviceType>();
  k_f_tab.modify_host(); k_f_tab.template sync<DeviceType>();
  k_df.modify_host(); k_df.template sync<DeviceType>();
  k_e2.modify_host(); k_e2.template sync<DeviceType>();
  k_f2.modify_host(); k_f2.template sync<DeviceType>();

  d_tabindex = k_tabindex.template view<DeviceType>();
  d_aat_k = k_aat_k.template view<DeviceType>();
  d_aat_theta0_1 = k_aat_theta0_1.template view<DeviceType>();
  d_aat_theta0_2 = k_aat_theta0_2.template view<DeviceType>();
  d_f_unspecified = k_f_unspecified.template view<DeviceType>();
  d_delta = k_delta.template view<DeviceType>();
  d_invdelta = k_invdelta.template view<DeviceType>();
  d_deltasq6 = k_deltasq6.template view<DeviceType>();
  d_e = k_e.template view<DeviceType>();
  d_de = k_de.template view<DeviceType>();
  d_f_tab = k_f_tab.template view<DeviceType>();
  d_df = k_df.template view<DeviceType>();
  d_e2 = k_e2.template view<DeviceType>();
  d_f2 = k_f2.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralTableCutKokkos<DeviceType>::init_style()
{
  DihedralTableCut::init_style();

  setup_tables();
}

/* ----------------------------------------------------------------------
   device version of DihedralTableCut::uf_lookup()
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralTableCutKokkos<DeviceType>::uf_lookup_kk(const int &type, const KK_FLOAT &x_in,
                                                   KK_FLOAT &u, KK_FLOAT &mdu) const
{
  const int tb = d_tabindex[type];
  const KK_FLOAT invdelta = d_invdelta[tb];
  const KK_FLOAT x_over_delta = x_in * invdelta;
  int i = static_cast<int>(x_over_delta);
  const KK_FLOAT b = x_over_delta - i;

  // the table is cyclic, so wrap both indices

  if (i >= tablength) i -= tablength;
  int ip1 = i + 1;
  if (ip1 >= tablength) ip1 -= tablength;

  if (tabstyle == LINEAR_STYLE) {

    // works even when the force column was not given in the table file

    u = d_e(tb,i) + b * d_de(tb,i);
    mdu = d_f_tab(tb,i) + b * d_df(tb,i);

  } else {

    const KK_FLOAT a = static_cast<KK_FLOAT>(1.0) - b;
    const KK_FLOAT deltasq6 = d_deltasq6[tb];

    u = a * d_e(tb,i) + b * d_e(tb,ip1) +
        ((a*a*a - a) * d_e2(tb,i) + (b*b*b - b) * d_e2(tb,ip1)) * deltasq6;

    if (d_f_unspecified[tb])

      // derivative of the energy spline, equation 3.3.5 of Numerical Recipes

      mdu = (d_e(tb,i) - d_e(tb,ip1)) * invdelta +
            ((static_cast<KK_FLOAT>(3.0)*a*a - static_cast<KK_FLOAT>(1.0)) * d_e2(tb,i) +
             (static_cast<KK_FLOAT>(1.0) - static_cast<KK_FLOAT>(3.0)*b*b) * d_e2(tb,ip1)) *
            d_delta[tb] / static_cast<KK_FLOAT>(6.0);
    else
      mdu = a * d_f_tab(tb,i) + b * d_f_tab(tb,ip1) +
            ((a*a*a - a) * d_f2(tb,i) + (b*b*b - b) * d_f2(tb,ip1)) * deltasq6;
  }
}


/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
          = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralTableCutKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                        KK_FLOAT &edihedral, KK_FLOAT *f1, KK_FLOAT *f3, KK_FLOAT *f4,
                        const KK_FLOAT &vb1x, const KK_FLOAT &vb1y, const KK_FLOAT &vb1z,
                        const KK_FLOAT &vb2x, const KK_FLOAT &vb2y, const KK_FLOAT &vb2z,
                        const KK_FLOAT &vb3x, const KK_FLOAT &vb3y, const KK_FLOAT &vb3z) const
{
  KK_FLOAT edihedralquarter;
  KK_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(edihedral);
      else {
        edihedralquarter = static_cast<KK_FLOAT>(0.25)*edihedral;
        if (i1 < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(edihedralquarter);
        if (i2 < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(edihedralquarter);
        if (i3 < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(edihedralquarter);
        if (i4 < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(edihedralquarter);
      }
    }
    if (eflag_atom) {
      edihedralquarter = static_cast<KK_FLOAT>(0.25)*edihedral;
      if (newton_bond || i1 < nlocal) v_eatom[i1] += static_cast<KK_ACC_FLOAT>(edihedralquarter);
      if (newton_bond || i2 < nlocal) v_eatom[i2] += static_cast<KK_ACC_FLOAT>(edihedralquarter);
      if (newton_bond || i3 < nlocal) v_eatom[i3] += static_cast<KK_ACC_FLOAT>(edihedralquarter);
      if (newton_bond || i4 < nlocal) v_eatom[i4] += static_cast<KK_ACC_FLOAT>(edihedralquarter);
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
        ev.v[0] += static_cast<KK_ACC_FLOAT>(v[0]);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(v[1]);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(v[2]);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(v[3]);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(v[4]);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(v[5]);
      } else {
        if (i1 < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
        }
        if (i2 < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
        }
        if (i3 < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
        }
        if (i4 < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i1 < nlocal) {
        v_vatom(i1,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
        v_vatom(i1,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
        v_vatom(i1,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
        v_vatom(i1,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
        v_vatom(i1,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
        v_vatom(i1,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
      }
      if (newton_bond || i2 < nlocal) {
        v_vatom(i2,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
        v_vatom(i2,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
        v_vatom(i2,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
        v_vatom(i2,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
        v_vatom(i2,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
        v_vatom(i2,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
      }
      if (newton_bond || i3 < nlocal) {
        v_vatom(i3,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
        v_vatom(i3,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
        v_vatom(i3,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
        v_vatom(i3,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
        v_vatom(i3,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
        v_vatom(i3,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
      }
      if (newton_bond || i4 < nlocal) {
        v_vatom(i4,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[0]);
        v_vatom(i4,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[1]);
        v_vatom(i4,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[2]);
        v_vatom(i4,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[3]);
        v_vatom(i4,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[4]);
        v_vatom(i4,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*v[5]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class DihedralTableCutKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class DihedralTableCutKokkos<LMPHostType>;
#endif
}

