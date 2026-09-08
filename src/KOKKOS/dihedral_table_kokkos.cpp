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

#include "dihedral_table_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "domain.h"
#include "math_const.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int LINEAR_STYLE = 0;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralTableKokkos<DeviceType>::DihedralTableKokkos(LAMMPS *lmp) : DihedralTable(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;


  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralTableKokkos<DeviceType>::~DihedralTableKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralTableKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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


  // copy domain data for the device-side minimum_image

  triclinic = domain->triclinic;
  xperiodic = domain->xperiodic;
  yperiodic = domain->yperiodic;
  zperiodic = domain->zperiodic;
  xprd = static_cast<KK_FLOAT>(domain->xprd);
  yprd = static_cast<KK_FLOAT>(domain->yprd);
  zprd = static_cast<KK_FLOAT>(domain->zprd);
  xprd_half = static_cast<KK_FLOAT>(domain->xprd_half);
  yprd_half = static_cast<KK_FLOAT>(domain->yprd_half);
  zprd_half = static_cast<KK_FLOAT>(domain->zprd_half);
  xy = static_cast<KK_FLOAT>(domain->xy);
  xz = static_cast<KK_FLOAT>(domain->xz);
  yz = static_cast<KK_FLOAT>(domain->yz);

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

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralTableCompute<1,1> >(0,ndihedrallist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralTableCompute<0,1> >(0,ndihedrallist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralTableCompute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralTableCompute<0,0> >(0,ndihedrallist),*this);
    }
  }

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
void DihedralTableKokkos<DeviceType>::operator()(TagDihedralTableCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {
  // The f array is atomic
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = dihedrallist(n,0);
  const int i2 = dihedrallist(n,1);
  const int i3 = dihedrallist(n,2);
  const int i4 = dihedrallist(n,3);
  const int type = dihedrallist(n,4);

  // ------ Step 1: compute the dihedral angle "phi" ------
  //
  // this is DihedralTable::Phi() inlined: it also leaves behind vb12,
  // vb23, vb34 and the two plane normals, which step 2 needs

  KK_FLOAT vb12[3],vb23[3],vb34[3],n123[3],n234[3];

  vb12[0] = x(i2,0) - x(i1,0);
  vb12[1] = x(i2,1) - x(i1,1);
  vb12[2] = x(i2,2) - x(i1,2);
  vb23[0] = x(i3,0) - x(i2,0);
  vb23[1] = x(i3,1) - x(i2,1);
  vb23[2] = x(i3,2) - x(i2,2);
  vb34[0] = x(i4,0) - x(i3,0);
  vb34[1] = x(i4,1) - x(i3,1);
  vb34[2] = x(i4,2) - x(i3,2);

  minimum_image(vb12[0],vb12[1],vb12[2]);
  minimum_image(vb23[0],vb23[1],vb23[2]);
  minimum_image(vb34[0],vb34[1],vb34[2]);

  MathSpecialKokkos::cross3(vb23,vb12,n123);
  MathSpecialKokkos::cross3(vb23,vb34,n234);

  KK_FLOAT n123_len = Kokkos::sqrt(MathSpecialKokkos::dot3(n123,n123));
  if (n123_len != static_cast<KK_FLOAT>(0.0)) {
    const KK_FLOAT s = static_cast<KK_FLOAT>(1.0)/n123_len;
    n123[0] *= s; n123[1] *= s; n123[2] *= s;
  }
  KK_FLOAT n234_len = Kokkos::sqrt(MathSpecialKokkos::dot3(n234,n234));
  if (n234_len != static_cast<KK_FLOAT>(0.0)) {
    const KK_FLOAT s = static_cast<KK_FLOAT>(1.0)/n234_len;
    n234[0] *= s; n234[1] *= s; n234[2] *= s;
  }

  KK_FLOAT cos_phi = -MathSpecialKokkos::dot3(n123,n234);
  if (cos_phi > static_cast<KK_FLOAT>(1.0)) cos_phi = static_cast<KK_FLOAT>(1.0);
  else if (cos_phi < static_cast<KK_FLOAT>(-1.0)) cos_phi = static_cast<KK_FLOAT>(-1.0);

  KK_FLOAT phi = Kokkos::acos(cos_phi);
  if (MathSpecialKokkos::dot3(n123,vb34) > static_cast<KK_FLOAT>(0.0)) {
    phi = -phi;
    phi += static_cast<KK_FLOAT>(MY_2PI);
  }

  // ------ Step 2: gradient of phi with respect to the atom positions ------

  KK_FLOAT dphi_dx1[3],dphi_dx2[3],dphi_dx3[3],dphi_dx4[3];
  KK_FLOAT proj12on23[3],proj34on23[3],perp12on23[3],perp34on23[3];

  const KK_FLOAT dot123 = MathSpecialKokkos::dot3(vb12,vb23);
  const KK_FLOAT dot234 = MathSpecialKokkos::dot3(vb23,vb34);
  const KK_FLOAT L23sqr = MathSpecialKokkos::dot3(vb23,vb23);
  const KK_FLOAT L23 = Kokkos::sqrt(L23sqr);

  KK_FLOAT inv_L23sqr = static_cast<KK_FLOAT>(0.0);
  KK_FLOAT inv_L23 = static_cast<KK_FLOAT>(0.0);
  if (L23sqr != static_cast<KK_FLOAT>(0.0)) {
    inv_L23sqr = static_cast<KK_FLOAT>(1.0)/L23sqr;
    inv_L23 = static_cast<KK_FLOAT>(1.0)/L23;
  }
  const KK_FLOAT neg_inv_L23 = -inv_L23;
  const KK_FLOAT dot123_over_L23sqr = dot123 * inv_L23sqr;
  const KK_FLOAT dot234_over_L23sqr = dot234 * inv_L23sqr;

  for (int d = 0; d < 3; ++d) {
    proj12on23[d] = vb23[d] * dot123_over_L23sqr;
    proj34on23[d] = vb23[d] * dot234_over_L23sqr;
    perp12on23[d] = vb12[d] - proj12on23[d];
    perp34on23[d] = vb34[d] - proj34on23[d];
  }

  const KK_FLOAT perp12on23_len = Kokkos::sqrt(MathSpecialKokkos::dot3(perp12on23,perp12on23));
  const KK_FLOAT perp34on23_len = Kokkos::sqrt(MathSpecialKokkos::dot3(perp34on23,perp34on23));

  KK_FLOAT inv_perp12on23 = static_cast<KK_FLOAT>(0.0);
  if (perp12on23_len != static_cast<KK_FLOAT>(0.0))
    inv_perp12on23 = static_cast<KK_FLOAT>(1.0)/perp12on23_len;
  KK_FLOAT inv_perp34on23 = static_cast<KK_FLOAT>(0.0);
  if (perp34on23_len != static_cast<KK_FLOAT>(0.0))
    inv_perp34on23 = static_cast<KK_FLOAT>(1.0)/perp34on23_len;

  for (int d = 0; d < 3; ++d) {
    dphi_dx1[d] = n123[d] * inv_perp12on23;
    dphi_dx4[d] = n234[d] * inv_perp34on23;
  }

  const KK_FLOAT proj12on23_len = dot123 * inv_L23;
  const KK_FLOAT proj34on23_len = dot234 * inv_L23;

  const KK_FLOAT dphi123_dx2_coef = neg_inv_L23 * (L23 + proj12on23_len);
  const KK_FLOAT dphi234_dx2_coef = inv_L23 * proj34on23_len;
  const KK_FLOAT dphi234_dx3_coef = neg_inv_L23 * (L23 + proj34on23_len);
  const KK_FLOAT dphi123_dx3_coef = inv_L23 * proj12on23_len;

  for (int d = 0; d < 3; ++d) {
    dphi_dx2[d] = dphi123_dx2_coef*dphi_dx1[d] + dphi234_dx2_coef*dphi_dx4[d];
    dphi_dx3[d] = dphi123_dx3_coef*dphi_dx1[d] + dphi234_dx3_coef*dphi_dx4[d];
  }

  // ------ Step 3: tabulated energy and force in the phi direction ------

  KK_FLOAT u = static_cast<KK_FLOAT>(0.0);
  KK_FLOAT m_du_dphi = static_cast<KK_FLOAT>(0.0);
  uf_lookup_kk(type,phi,u,m_du_dphi);

  KK_FLOAT edihedral = 0;
  if (eflag) edihedral = u;

  // ------ Step 4: chain rule to get the force in real space ------

  KK_FLOAT f1[3],f2[3],f3[3],f4[3];
  for (int d = 0; d < 3; ++d) {
    f1[d] = m_du_dphi * dphi_dx1[d];
    f2[d] = m_du_dphi * dphi_dx2[d];
    f3[d] = m_du_dphi * dphi_dx3[d];
    f4[d] = m_du_dphi * dphi_dx4[d];
  }

  // apply force to each of 4 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += static_cast<KK_ACC_FLOAT>(f1[0]);
    a_f(i1,1) += static_cast<KK_ACC_FLOAT>(f1[1]);
    a_f(i1,2) += static_cast<KK_ACC_FLOAT>(f1[2]);
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) += static_cast<KK_ACC_FLOAT>(f2[0]);
    a_f(i2,1) += static_cast<KK_ACC_FLOAT>(f2[1]);
    a_f(i2,2) += static_cast<KK_ACC_FLOAT>(f2[2]);
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += static_cast<KK_ACC_FLOAT>(f3[0]);
    a_f(i3,1) += static_cast<KK_ACC_FLOAT>(f3[1]);
    a_f(i3,2) += static_cast<KK_ACC_FLOAT>(f3[2]);
  }

  if (NEWTON_BOND || i4 < nlocal) {
    a_f(i4,0) += static_cast<KK_ACC_FLOAT>(f4[0]);
    a_f(i4,1) += static_cast<KK_ACC_FLOAT>(f4[1]);
    a_f(i4,2) += static_cast<KK_ACC_FLOAT>(f4[2]);
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,edihedral,f1,f3,f4,
             -vb12[0],-vb12[1],-vb12[2],
             vb23[0],vb23[1],vb23[2],
             vb34[0],vb34[1],vb34[2]);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralTableKokkos<DeviceType>::operator()(TagDihedralTableCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralTableCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ----------------------------------------------------------------------
   copy the tabulated data built on the host by compute_table() into
   device views.  this runs once per setup, after every coeff() call, so
   restarts and repeated dihedral_coeff commands are covered
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralTableKokkos<DeviceType>::setup_tables()
{
  const int n = atom->ndihedraltypes;

  k_tabindex = DAT::tdual_int_1d("DihedralTable::tabindex",n+1);
  k_f_unspecified = DAT::tdual_int_1d("DihedralTable::f_unspecified",ntables);
  k_delta = DAT::tdual_kkfloat_1d("DihedralTable::delta",ntables);
  k_invdelta = DAT::tdual_kkfloat_1d("DihedralTable::invdelta",ntables);
  k_deltasq6 = DAT::tdual_kkfloat_1d("DihedralTable::deltasq6",ntables);

  k_e = DAT::tdual_kkfloat_2d("DihedralTable::e",ntables,tablength);
  k_de = DAT::tdual_kkfloat_2d("DihedralTable::de",ntables,tablength);
  k_f_tab = DAT::tdual_kkfloat_2d("DihedralTable::f",ntables,tablength);
  k_df = DAT::tdual_kkfloat_2d("DihedralTable::df",ntables,tablength);
  k_e2 = DAT::tdual_kkfloat_2d("DihedralTable::e2",ntables,tablength);
  k_f2 = DAT::tdual_kkfloat_2d("DihedralTable::f2",ntables,tablength);

  for (int i = 1; i <= n; i++) k_tabindex.view_host()(i) = tabindex[i];

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
void DihedralTableKokkos<DeviceType>::init_style()
{
  DihedralTable::init_style();

  setup_tables();
}

/* ----------------------------------------------------------------------
   device version of DihedralTable::uf_lookup()
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralTableKokkos<DeviceType>::uf_lookup_kk(const int &type, const KK_FLOAT &x_in,
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
   device version of Domain::minimum_image()
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralTableKokkos<DeviceType>::minimum_image(KK_FLOAT &dx, KK_FLOAT &dy, KK_FLOAT &dz) const
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (Kokkos::abs(dx) > xprd_half) {
        if (dx < static_cast<KK_FLOAT>(0.0)) dx += xprd;
        else dx -= xprd;
      }
    }
    if (yperiodic) {
      if (Kokkos::abs(dy) > yprd_half) {
        if (dy < static_cast<KK_FLOAT>(0.0)) dy += yprd;
        else dy -= yprd;
      }
    }
    if (zperiodic) {
      if (Kokkos::abs(dz) > zprd_half) {
        if (dz < static_cast<KK_FLOAT>(0.0)) dz += zprd;
        else dz -= zprd;
      }
    }
  } else {
    if (zperiodic) {
      if (Kokkos::abs(dz) > zprd_half) {
        if (dz < static_cast<KK_FLOAT>(0.0)) {
          dz += zprd; dy += yz; dx += xz;
        } else {
          dz -= zprd; dy -= yz; dx -= xz;
        }
      }
    }
    if (yperiodic) {
      if (Kokkos::abs(dy) > yprd_half) {
        if (dy < static_cast<KK_FLOAT>(0.0)) {
          dy += yprd; dx += xy;
        } else {
          dy -= yprd; dx -= xy;
        }
      }
    }
    if (xperiodic) {
      if (Kokkos::abs(dx) > xprd_half) {
        if (dx < static_cast<KK_FLOAT>(0.0)) dx += xprd;
        else dx -= xprd;
      }
    }
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
void DihedralTableKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
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
template class DihedralTableKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class DihedralTableKokkos<LMPHostType>;
#endif
}

