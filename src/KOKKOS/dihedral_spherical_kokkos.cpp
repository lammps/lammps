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
   Contributing author: Axel Kohlmeyer (Temple U)
   [ based on dihedral_nharmonic_kokkos.cpp and dihedral_spherical.cpp ]
------------------------------------------------------------------------- */

#include "dihedral_spherical_kokkos.h"

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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralSphericalKokkos<DeviceType>::DihedralSphericalKokkos(LAMMPS *lmp) : DihedralSpherical(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = DAT::tdual_int_scalar("Dihedral:warning_flag");
  d_warning_flag = k_warning_flag.view<DeviceType>();
  h_warning_flag = k_warning_flag.view_host();

  centroidstressflag = CENTROID_NOTAVAIL;

  allocated_kokkos = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralSphericalKokkos<DeviceType>::~DihedralSphericalKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralSphericalKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    if ((int)k_eatom.extent(0) < maxeatom) {
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"dihedral:eatom");
      d_eatom = k_eatom.view<DeviceType>();
    } else Kokkos::deep_copy(d_eatom,0.0);
  }
  if (vflag_atom) {
    if ((int)k_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"dihedral:vatom");
      d_vatom = k_vatom.view<DeviceType>();
    } else Kokkos::deep_copy(d_vatom,0.0);
  }

  k_nterms.template sync<DeviceType>();
  k_Ccoeff.template sync<DeviceType>();
  k_phi_mult.template sync<DeviceType>();
  k_phi_shift.template sync<DeviceType>();
  k_phi_offset.template sync<DeviceType>();
  k_theta1_mult.template sync<DeviceType>();
  k_theta1_shift.template sync<DeviceType>();
  k_theta1_offset.template sync<DeviceType>();
  k_theta2_mult.template sync<DeviceType>();
  k_theta2_shift.template sync<DeviceType>();
  k_theta2_offset.template sync<DeviceType>();

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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralSphericalCompute<1,1> >(0,ndihedrallist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralSphericalCompute<0,1> >(0,ndihedrallist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralSphericalCompute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralSphericalCompute<0,0> >(0,ndihedrallist),*this);
    }
  }

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.sync_host();
  if (h_warning_flag())
    error->warning(FLERR,"Dihedral problem");

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
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
KK_FLOAT DihedralSphericalKokkos<DeviceType>::CalcGeneralizedForcesKK(
    int type, KK_FLOAT phi, KK_FLOAT theta1, KK_FLOAT theta2,
    KK_FLOAT &m_du_dth1, KK_FLOAT &m_du_dth2, KK_FLOAT &m_du_dphi) const
{
  KK_FLOAT energy = 0.0;
  m_du_dphi = 0.0;
  m_du_dth1 = 0.0;
  m_du_dth2 = 0.0;

  const int nt = d_nterms[type];

  for (int j = 0; j < nt; j++) {
    KK_FLOAT cp = 1.0;
    KK_FLOAT sp = 0.0;
    const KK_FLOAT pm = d_phi_mult(type,j);
    if (pm != 0.0) {
      const KK_FLOAT p = pm * (phi - d_phi_shift(type,j));
      cp = cos(p);
      sp = sin(p);
    }

    KK_FLOAT ct1 = 1.0;
    KK_FLOAT st1 = 0.0;
    const KK_FLOAT t1m = d_theta1_mult(type,j);
    if (t1m != 0.0) {
      const KK_FLOAT t1 = t1m * (theta1 - d_theta1_shift(type,j));
      ct1 = cos(t1);
      st1 = sin(t1);
    }

    KK_FLOAT ct2 = 1.0;
    KK_FLOAT st2 = 0.0;
    const KK_FLOAT t2m = d_theta2_mult(type,j);
    if (t2m != 0.0) {
      const KK_FLOAT t2 = t2m * (theta2 - d_theta2_shift(type,j));
      ct2 = cos(t2);
      st2 = sin(t2);
    }

    const KK_FLOAT C = d_Ccoeff(type,j);
    const KK_FLOAT poff = d_phi_offset(type,j);
    const KK_FLOAT t1off = d_theta1_offset(type,j);
    const KK_FLOAT t2off = d_theta2_offset(type,j);

    energy += C * (poff - cp) * (t1off - ct1) * (t2off - ct2);

    m_du_dphi += -C * sp * pm * (t1off - ct1) * (t2off - ct2);
    m_du_dth1 += -C * (poff - cp) * st1 * t1m * (t2off - ct2);
    m_du_dth2 += -C * (poff - cp) * (t1off - ct1) * st2 * t2m;
  }

  return energy;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralSphericalKokkos<DeviceType>::operator()(TagDihedralSphericalCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = dihedrallist(n,0);
  const int i2 = dihedrallist(n,1);
  const int i3 = dihedrallist(n,2);
  const int i4 = dihedrallist(n,3);
  const int type = dihedrallist(n,4);

  static constexpr int g_dim = 3;

  // Bond vectors
  KK_FLOAT vb12[g_dim], vb23[g_dim], vb34[g_dim];
  for (int d = 0; d < g_dim; ++d) {
    vb12[d] = x(i2,d) - x(i1,d);
    vb23[d] = x(i3,d) - x(i2,d);
    vb34[d] = x(i4,d) - x(i3,d);
  }

  // Compute normals to planes 1,2,3 and 2,3,4
  // n123 = vb23 x vb12
  // n234 = vb23 x vb34
  KK_FLOAT n123[g_dim], n234[g_dim];
  n123[0] = vb23[1]*vb12[2] - vb23[2]*vb12[1];
  n123[1] = vb23[2]*vb12[0] - vb23[0]*vb12[2];
  n123[2] = vb23[0]*vb12[1] - vb23[1]*vb12[0];

  n234[0] = vb23[1]*vb34[2] - vb23[2]*vb34[1];
  n234[1] = vb23[2]*vb34[0] - vb23[0]*vb34[2];
  n234[2] = vb23[0]*vb34[1] - vb23[1]*vb34[0];

  // Normalize n123 and n234
  KK_FLOAT inv_scale;
  inv_scale = sqrt(n123[0]*n123[0] + n123[1]*n123[1] + n123[2]*n123[2]);
  if (inv_scale > 0.0) {
    KK_FLOAT scale = 1.0/inv_scale;
    n123[0] *= scale; n123[1] *= scale; n123[2] *= scale;
  }
  inv_scale = sqrt(n234[0]*n234[0] + n234[1]*n234[1] + n234[2]*n234[2]);
  if (inv_scale > 0.0) {
    KK_FLOAT scale = 1.0/inv_scale;
    n234[0] *= scale; n234[1] *= scale; n234[2] *= scale;
  }

  // Dihedral angle phi
  KK_FLOAT cos_phi = -(n123[0]*n234[0] + n123[1]*n234[1] + n123[2]*n234[2]);
  if (cos_phi > 1.0) cos_phi = 1.0;
  else if (cos_phi < -1.0) cos_phi = -1.0;
  KK_FLOAT phi = acos(cos_phi);

  // Determine sign: if n123 . vb34 > 0 => negative dihedral
  KK_FLOAT n123_dot_vb34 = n123[0]*vb34[0] + n123[1]*vb34[1] + n123[2]*vb34[2];
  if (n123_dot_vb34 > 0.0) {
    phi = -phi;
    phi += MY_2PI;
  }

  // Dot products needed for bond lengths
  KK_FLOAT dot123 = vb12[0]*vb23[0] + vb12[1]*vb23[1] + vb12[2]*vb23[2];
  KK_FLOAT dot234 = vb23[0]*vb34[0] + vb23[1]*vb34[1] + vb23[2]*vb34[2];

  KK_FLOAT L12sqr = vb12[0]*vb12[0] + vb12[1]*vb12[1] + vb12[2]*vb12[2];
  KK_FLOAT L23sqr = vb23[0]*vb23[0] + vb23[1]*vb23[1] + vb23[2]*vb23[2];
  KK_FLOAT L34sqr = vb34[0]*vb34[0] + vb34[1]*vb34[1] + vb34[2]*vb34[2];

  KK_FLOAT L12 = sqrt(L12sqr);
  KK_FLOAT L23 = sqrt(L23sqr);
  KK_FLOAT L34 = sqrt(L34sqr);

  KK_FLOAT inv_L12sqr = (L12sqr != 0.0) ? 1.0/L12sqr : 0.0;
  KK_FLOAT inv_L12    = (L12sqr != 0.0) ? 1.0/L12    : 0.0;
  KK_FLOAT inv_L23sqr = (L23sqr != 0.0) ? 1.0/L23sqr : 0.0;
  KK_FLOAT inv_L23    = (L23sqr != 0.0) ? 1.0/L23    : 0.0;
  KK_FLOAT inv_L34sqr = (L34sqr != 0.0) ? 1.0/L34sqr : 0.0;
  KK_FLOAT inv_L34    = (L34sqr != 0.0) ? 1.0/L34    : 0.0;

  KK_FLOAT neg_inv_L23 = -inv_L23;

  KK_FLOAT dot123_over_L23sqr = dot123 * inv_L23sqr;
  KK_FLOAT dot234_over_L23sqr = dot234 * inv_L23sqr;
  KK_FLOAT dot123_over_L12sqr = dot123 * inv_L12sqr;
  KK_FLOAT dot234_over_L34sqr = dot234 * inv_L34sqr;

  // Projection and perpendicular vectors
  KK_FLOAT proj12on23[g_dim], perp12on23[g_dim];
  KK_FLOAT proj34on23[g_dim], perp34on23[g_dim];
  KK_FLOAT proj23on12[g_dim], perp23on12[g_dim];
  KK_FLOAT proj23on34[g_dim], perp23on34[g_dim];

  for (int d = 0; d < g_dim; ++d) {
    proj12on23[d] = vb23[d] * dot123_over_L23sqr;
    proj34on23[d] = vb23[d] * dot234_over_L23sqr;
    perp12on23[d] = vb12[d] - proj12on23[d];
    perp34on23[d] = vb34[d] - proj34on23[d];

    proj23on12[d] = vb12[d] * dot123_over_L12sqr;
    proj23on34[d] = vb34[d] * dot234_over_L34sqr;
    perp23on12[d] = vb23[d] - proj23on12[d];
    perp23on34[d] = vb23[d] - proj23on34[d];
  }

  KK_FLOAT perp12on23_len = sqrt(perp12on23[0]*perp12on23[0] + perp12on23[1]*perp12on23[1] + perp12on23[2]*perp12on23[2]);
  KK_FLOAT perp34on23_len = sqrt(perp34on23[0]*perp34on23[0] + perp34on23[1]*perp34on23[1] + perp34on23[2]*perp34on23[2]);
  KK_FLOAT perp23on12_len = sqrt(perp23on12[0]*perp23on12[0] + perp23on12[1]*perp23on12[1] + perp23on12[2]*perp23on12[2]);
  KK_FLOAT perp23on34_len = sqrt(perp23on34[0]*perp23on34[0] + perp23on34[1]*perp23on34[1] + perp23on34[2]*perp23on34[2]);

  KK_FLOAT inv_perp12on23 = (perp12on23_len != 0.0) ? 1.0/perp12on23_len : 0.0;
  KK_FLOAT inv_perp34on23 = (perp34on23_len != 0.0) ? 1.0/perp34on23_len : 0.0;
  KK_FLOAT inv_perp23on12 = (perp23on12_len != 0.0) ? 1.0/perp23on12_len : 0.0;
  KK_FLOAT inv_perp23on34 = (perp23on34_len != 0.0) ? 1.0/perp23on34_len : 0.0;

  // Gradients of phi
  KK_FLOAT dphi_dx1[g_dim], dphi_dx2[g_dim], dphi_dx3[g_dim], dphi_dx4[g_dim];

  for (int d = 0; d < g_dim; ++d) {
    dphi_dx1[d] = n123[d] * inv_perp12on23;
    dphi_dx4[d] = n234[d] * inv_perp34on23;
  }

  KK_FLOAT proj12on23_len = dot123 * inv_L23;
  KK_FLOAT proj34on23_len = dot234 * inv_L23;

  KK_FLOAT dphi123_dx2_coef = neg_inv_L23 * (L23 + proj12on23_len);
  KK_FLOAT dphi234_dx2_coef = inv_L23 * proj34on23_len;
  KK_FLOAT dphi234_dx3_coef = neg_inv_L23 * (L23 + proj34on23_len);
  KK_FLOAT dphi123_dx3_coef = inv_L23 * proj12on23_len;

  for (int d = 0; d < g_dim; ++d) {
    dphi_dx2[d] = dphi123_dx2_coef * dphi_dx1[d] + dphi234_dx2_coef * dphi_dx4[d];
    dphi_dx3[d] = dphi123_dx3_coef * dphi_dx1[d] + dphi234_dx3_coef * dphi_dx4[d];
  }

  // Gradients of theta1 and theta2
  KK_FLOAT dth1_dx1[g_dim], dth1_dx2[g_dim], dth1_dx3[g_dim];
  KK_FLOAT dth2_dx2[g_dim], dth2_dx3[g_dim], dth2_dx4[g_dim];

  KK_FLOAT coeff_dth1_dx1 = -inv_perp23on12 * inv_L12;
  KK_FLOAT coeff_dth1_dx3 = inv_perp12on23 * inv_L23;
  KK_FLOAT coeff_dth2_dx2 = -inv_perp34on23 * inv_L23;
  KK_FLOAT coeff_dth2_dx4 = inv_perp23on34 * inv_L34;

  for (int d = 0; d < g_dim; ++d) {
    dth1_dx1[d] = perp23on12[d] * coeff_dth1_dx1;
    dth1_dx3[d] = perp12on23[d] * coeff_dth1_dx3;
    dth1_dx2[d] = -(dth1_dx1[d] + dth1_dx3[d]);

    dth2_dx2[d] = perp34on23[d] * coeff_dth2_dx2;
    dth2_dx4[d] = perp23on34[d] * coeff_dth2_dx4;
    dth2_dx3[d] = -(dth2_dx2[d] + dth2_dx4[d]);
  }

  // Bond angles theta1 and theta2
  KK_FLOAT ct1 = -dot123 * inv_L12 * inv_L23;
  if (ct1 < -1.0) ct1 = -1.0;
  else if (ct1 > 1.0) ct1 = 1.0;
  KK_FLOAT theta1 = acos(ct1);

  KK_FLOAT ct2 = -dot234 * inv_L23 * inv_L34;
  if (ct2 < -1.0) ct2 = -1.0;
  else if (ct2 > 1.0) ct2 = 1.0;
  KK_FLOAT theta2 = acos(ct2);

  // Generalized forces
  KK_FLOAT m_du_dth1 = 0.0, m_du_dth2 = 0.0, m_du_dphi = 0.0;
  KK_FLOAT edihedral = CalcGeneralizedForcesKK(type, phi, theta1, theta2,
                                               m_du_dth1, m_du_dth2, m_du_dphi);

  // Forces on all 4 atoms
  KK_FLOAT f1[g_dim], f2[g_dim], f3[g_dim], f4[g_dim];
  for (int d = 0; d < g_dim; ++d) {
    f1[d] = m_du_dphi * dphi_dx1[d] + m_du_dth1 * dth1_dx1[d];
    f2[d] = m_du_dphi * dphi_dx2[d] + m_du_dth1 * dth1_dx2[d] + m_du_dth2 * dth2_dx2[d];
    f3[d] = m_du_dphi * dphi_dx3[d] + m_du_dth1 * dth1_dx3[d] + m_du_dth2 * dth2_dx3[d];
    f4[d] = m_du_dphi * dphi_dx4[d] + m_du_dth2 * dth2_dx4[d];
  }

  // Apply force to each of 4 atoms

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

  // ev_tally uses vb12 as vb1, vb23 as vb2, vb34 as vb3
  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,edihedral,f1,f3,f4,
             vb12[0],vb12[1],vb12[2],vb23[0],vb23[1],vb23[2],vb34[0],vb34[1],vb34[2]);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralSphericalKokkos<DeviceType>::operator()(TagDihedralSphericalCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralSphericalCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralSphericalKokkos<DeviceType>::allocate_kokkos()
{
  int n = atom->ndihedraltypes;

  if (!allocated_kokkos) {
    k_nterms = DAT::tdual_int_1d("DihedralSpherical::nterms",n+1);
    k_Ccoeff = DAT::tdual_kkfloat_2d("DihedralSpherical::Ccoeff",n+1,nterms_max);
    k_phi_mult = DAT::tdual_kkfloat_2d("DihedralSpherical::phi_mult",n+1,nterms_max);
    k_phi_shift = DAT::tdual_kkfloat_2d("DihedralSpherical::phi_shift",n+1,nterms_max);
    k_phi_offset = DAT::tdual_kkfloat_2d("DihedralSpherical::phi_offset",n+1,nterms_max);
    k_theta1_mult = DAT::tdual_kkfloat_2d("DihedralSpherical::theta1_mult",n+1,nterms_max);
    k_theta1_shift = DAT::tdual_kkfloat_2d("DihedralSpherical::theta1_shift",n+1,nterms_max);
    k_theta1_offset = DAT::tdual_kkfloat_2d("DihedralSpherical::theta1_offset",n+1,nterms_max);
    k_theta2_mult = DAT::tdual_kkfloat_2d("DihedralSpherical::theta2_mult",n+1,nterms_max);
    k_theta2_shift = DAT::tdual_kkfloat_2d("DihedralSpherical::theta2_shift",n+1,nterms_max);
    k_theta2_offset = DAT::tdual_kkfloat_2d("DihedralSpherical::theta2_offset",n+1,nterms_max);
  } else {
    k_nterms.resize(n+1);
    k_Ccoeff.resize(n+1,nterms_max);
    k_phi_mult.resize(n+1,nterms_max);
    k_phi_shift.resize(n+1,nterms_max);
    k_phi_offset.resize(n+1,nterms_max);
    k_theta1_mult.resize(n+1,nterms_max);
    k_theta1_shift.resize(n+1,nterms_max);
    k_theta1_offset.resize(n+1,nterms_max);
    k_theta2_mult.resize(n+1,nterms_max);
    k_theta2_shift.resize(n+1,nterms_max);
    k_theta2_offset.resize(n+1,nterms_max);
  }

  d_nterms = k_nterms.template view<DeviceType>();
  d_Ccoeff = k_Ccoeff.template view<DeviceType>();
  d_phi_mult = k_phi_mult.template view<DeviceType>();
  d_phi_shift = k_phi_shift.template view<DeviceType>();
  d_phi_offset = k_phi_offset.template view<DeviceType>();
  d_theta1_mult = k_theta1_mult.template view<DeviceType>();
  d_theta1_shift = k_theta1_shift.template view<DeviceType>();
  d_theta1_offset = k_theta1_offset.template view<DeviceType>();
  d_theta2_mult = k_theta2_mult.template view<DeviceType>();
  d_theta2_shift = k_theta2_shift.template view<DeviceType>();
  d_theta2_offset = k_theta2_offset.template view<DeviceType>();

  allocated_kokkos = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralSphericalKokkos<DeviceType>::allocate()
{
  DihedralSpherical::allocate();
  // Kokkos views will be allocated in coeff() once nterms_max is known
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralSphericalKokkos<DeviceType>::coeff(int narg, char **arg)
{
  DihedralSpherical::coeff(narg, arg);
  allocate_kokkos();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->ndihedraltypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_nterms.view_host()[i] = nterms[i];
    for (int j = 0; j < nterms[i]; j++) {
      k_Ccoeff.view_host()(i,j) = Ccoeff[i][j];
      k_phi_mult.view_host()(i,j) = phi_mult[i][j];
      k_phi_shift.view_host()(i,j) = phi_shift[i][j];
      k_phi_offset.view_host()(i,j) = phi_offset[i][j];
      k_theta1_mult.view_host()(i,j) = theta1_mult[i][j];
      k_theta1_shift.view_host()(i,j) = theta1_shift[i][j];
      k_theta1_offset.view_host()(i,j) = theta1_offset[i][j];
      k_theta2_mult.view_host()(i,j) = theta2_mult[i][j];
      k_theta2_shift.view_host()(i,j) = theta2_shift[i][j];
      k_theta2_offset.view_host()(i,j) = theta2_offset[i][j];
    }
  }

  k_nterms.modify_host();
  k_Ccoeff.modify_host();
  k_phi_mult.modify_host();
  k_phi_shift.modify_host();
  k_phi_offset.modify_host();
  k_theta1_mult.modify_host();
  k_theta1_shift.modify_host();
  k_theta1_offset.modify_host();
  k_theta2_mult.modify_host();
  k_theta2_shift.modify_host();
  k_theta2_offset.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralSphericalKokkos<DeviceType>::read_restart(FILE *fp)
{
  DihedralSpherical::read_restart(fp);
  allocate_kokkos();

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_nterms.view_host()[i] = nterms[i];
    for (int j = 0; j < nterms[i]; j++) {
      k_Ccoeff.view_host()(i,j) = Ccoeff[i][j];
      k_phi_mult.view_host()(i,j) = phi_mult[i][j];
      k_phi_shift.view_host()(i,j) = phi_shift[i][j];
      k_phi_offset.view_host()(i,j) = phi_offset[i][j];
      k_theta1_mult.view_host()(i,j) = theta1_mult[i][j];
      k_theta1_shift.view_host()(i,j) = theta1_shift[i][j];
      k_theta1_offset.view_host()(i,j) = theta1_offset[i][j];
      k_theta2_mult.view_host()(i,j) = theta2_mult[i][j];
      k_theta2_shift.view_host()(i,j) = theta2_shift[i][j];
      k_theta2_offset.view_host()(i,j) = theta2_offset[i][j];
    }
  }

  k_nterms.modify_host();
  k_Ccoeff.modify_host();
  k_phi_mult.modify_host();
  k_phi_shift.modify_host();
  k_phi_offset.modify_host();
  k_theta1_mult.modify_host();
  k_theta1_shift.modify_host();
  k_theta1_offset.modify_host();
  k_theta2_mult.modify_host();
  k_theta2_shift.modify_host();
  k_theta2_offset.modify_host();
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
void DihedralSphericalKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
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
template class DihedralSphericalKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class DihedralSphericalKokkos<LMPHostType>;
#endif
}
