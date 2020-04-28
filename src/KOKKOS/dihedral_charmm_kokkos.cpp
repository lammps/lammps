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

#include "dihedral_charmm_kokkos.h"
#include <cmath>
#include <cstdlib>
#include "atom_kokkos.h"
#include "comm.h"
#include "neighbor_kokkos.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralCharmmKokkos<DeviceType>::DihedralCharmmKokkos(LAMMPS *lmp) : DihedralCharmm(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK | TYPE_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = Kokkos::DualView<int,DeviceType>("Dihedral:warning_flag");
  d_warning_flag = k_warning_flag.template view<DeviceType>();
  h_warning_flag = k_warning_flag.h_view;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralCharmmKokkos<DeviceType>::~DihedralCharmmKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCharmmKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // insure pair->ev_tally() will use 1-4 virial contribution

  if (weightflag && vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    //if(k_eatom.extent(0)<maxeatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"dihedral:eatom");
      d_eatom = k_eatom.template view<KKDeviceType>();
      k_eatom_pair = Kokkos::DualView<E_FLOAT*,Kokkos::LayoutRight,KKDeviceType>("dihedral:eatom_pair",maxeatom);
      d_eatom_pair = k_eatom.template view<KKDeviceType>();
    //}
  }
  if (vflag_atom) {
    //if(k_vatom.extent(0)<maxvatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"dihedral:vatom");
      d_vatom = k_vatom.template view<KKDeviceType>();
      k_vatom_pair = Kokkos::DualView<F_FLOAT*[6],Kokkos::LayoutRight,KKDeviceType>("dihedral:vatom_pair",maxvatom);
      d_vatom_pair = k_vatom.template view<KKDeviceType>();
    //}
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  atomtype = atomKK->k_type.view<DeviceType>();
  neighborKK->k_dihedrallist.template sync<DeviceType>();
  dihedrallist = neighborKK->k_dihedrallist.view<DeviceType>();
  int ndihedrallist = neighborKK->ndihedrallist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;
  qqrd2e = force->qqrd2e;

  h_warning_flag() = 0;
  k_warning_flag.template modify<LMPHostType>();
  k_warning_flag.template sync<DeviceType>();

  copymode = 1;

  // loop over neighbors of my atoms

  EVM_FLOAT evm;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralCharmmCompute<1,1> >(0,ndihedrallist),*this,evm);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralCharmmCompute<0,1> >(0,ndihedrallist),*this,evm);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralCharmmCompute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralCharmmCompute<0,0> >(0,ndihedrallist),*this);
    }
  }

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.template sync<LMPHostType>();
  if (h_warning_flag())
    error->warning(FLERR,"Dihedral problem",0);

  if (eflag_global) {
    energy += evm.emol;
    force->pair->eng_vdwl += evm.evdwl;
    force->pair->eng_coul += evm.ecoul;
  }
  if (vflag_global) {
    virial[0] += evm.v[0];
    virial[1] += evm.v[1];
    virial[2] += evm.v[2];
    virial[3] += evm.v[3];
    virial[4] += evm.v[4];
    virial[5] += evm.v[5];

    force->pair->virial[0] += evm.vp[0];
    force->pair->virial[1] += evm.vp[1];
    force->pair->virial[2] += evm.vp[2];
    force->pair->virial[3] += evm.vp[3];
    force->pair->virial[4] += evm.vp[4];
    force->pair->virial[5] += evm.vp[5];
  }

  // don't yet have dualviews for eatom and vatom in pair_kokkos,
  //  so need to manually copy these to pair style

  int n = nlocal;
  if (newton_bond) n += atom->nghost;

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();

    k_eatom_pair.template modify<DeviceType>();
    k_eatom_pair.template sync<LMPHostType>();
    for (int i = 0; i < n; i++)
      force->pair->eatom[i] += k_eatom_pair.h_view(i);
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();

    k_vatom_pair.template modify<DeviceType>();
    k_vatom_pair.template sync<LMPHostType>();
    for (int i = 0; i < n; i++) {
      force->pair->vatom[i][0] += k_vatom_pair.h_view(i,0);
      force->pair->vatom[i][1] += k_vatom_pair.h_view(i,1);
      force->pair->vatom[i][2] += k_vatom_pair.h_view(i,2);
      force->pair->vatom[i][3] += k_vatom_pair.h_view(i,3);
      force->pair->vatom[i][4] += k_vatom_pair.h_view(i,4);
      force->pair->vatom[i][5] += k_vatom_pair.h_view(i,5);
    }
  }

  copymode = 0;
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void DihedralCharmmKokkos<DeviceType>::operator()(TagDihedralCharmmCompute<NEWTON_BOND,EVFLAG>, const int &n, EVM_FLOAT& evm) const {

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

  const F_FLOAT ax = vb1y*vb2zm - vb1z*vb2ym;
  const F_FLOAT ay = vb1z*vb2xm - vb1x*vb2zm;
  const F_FLOAT az = vb1x*vb2ym - vb1y*vb2xm;
  const F_FLOAT bx = vb3y*vb2zm - vb3z*vb2ym;
  const F_FLOAT by = vb3z*vb2xm - vb3x*vb2zm;
  const F_FLOAT bz = vb3x*vb2ym - vb3y*vb2xm;

  const F_FLOAT rasq = ax*ax + ay*ay + az*az;
  const F_FLOAT rbsq = bx*bx + by*by + bz*bz;
  const F_FLOAT rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
  const F_FLOAT rg = sqrt(rgsq);

  F_FLOAT rginv,ra2inv,rb2inv;
  rginv = ra2inv = rb2inv = 0.0;
  if (rg > 0) rginv = 1.0/rg;
  if (rasq > 0) ra2inv = 1.0/rasq;
  if (rbsq > 0) rb2inv = 1.0/rbsq;
  const F_FLOAT rabinv = sqrt(ra2inv*rb2inv);

  F_FLOAT c = (ax*bx + ay*by + az*bz)*rabinv;
  F_FLOAT s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

    // error check

  if ((c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) && !d_warning_flag())
    Kokkos::atomic_fetch_add(&d_warning_flag(),1);

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  const int m = d_multiplicity[type];
  F_FLOAT p = 1.0;
  F_FLOAT ddf1,df1;
  ddf1 = df1 = 0.0;

  for (int i = 0; i < m; i++) {
    ddf1 = p*c - df1*s;
    df1 = p*s + df1*c;
    p = ddf1;
  }

  p = p*d_cos_shift[type] + df1*d_sin_shift[type];
  df1 = df1*d_cos_shift[type] - ddf1*d_sin_shift[type];
  df1 *= -m;
  p += 1.0;

  if (m == 0) {
    p = 1.0 + d_cos_shift[type];
    df1 = 0.0;
  }

  E_FLOAT edihedral = 0.0;
  if (eflag) edihedral = d_k[type] * p;

  const F_FLOAT fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
  const F_FLOAT hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
  const F_FLOAT fga = fg*ra2inv*rginv;
  const F_FLOAT hgb = hg*rb2inv*rginv;
  const F_FLOAT gaa = -ra2inv*rg;
  const F_FLOAT gbb = rb2inv*rg;

  const F_FLOAT dtfx = gaa*ax;
  const F_FLOAT dtfy = gaa*ay;
  const F_FLOAT dtfz = gaa*az;
  const F_FLOAT dtgx = fga*ax - hgb*bx;
  const F_FLOAT dtgy = fga*ay - hgb*by;
  const F_FLOAT dtgz = fga*az - hgb*bz;
  const F_FLOAT dthx = gbb*bx;
  const F_FLOAT dthy = gbb*by;
  const F_FLOAT dthz = gbb*bz;

  const F_FLOAT df = -d_k[type] * df1;

  const F_FLOAT sx2 = df*dtgx;
  const F_FLOAT sy2 = df*dtgy;
  const F_FLOAT sz2 = df*dtgz;

  F_FLOAT f1[3],f2[3],f3[3],f4[3];
  f1[0] = df*dtfx;
  f1[1] = df*dtfy;
  f1[2] = df*dtfz;

  f2[0] = sx2 - f1[0];
  f2[1] = sy2 - f1[1];
  f2[2] = sz2 - f1[2];

  f4[0] = df*dthx;
  f4[1] = df*dthy;
  f4[2] = df*dthz;

  f3[0] = -sx2 - f4[0];
  f3[1] = -sy2 - f4[1];
  f3[2] = -sz2 - f4[2];

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
    ev_tally(evm,i1,i2,i3,i4,edihedral,f1,f3,f4,
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);

  // 1-4 LJ and Coulomb interactions
  // tally energy/virial in pair, using newton_bond as newton flag

  if (d_weight[type] > 0.0) {
    const int itype = atomtype[i1];
    const int jtype = atomtype[i4];

    const F_FLOAT delx = x(i1,0) - x(i4,0);
    const F_FLOAT dely = x(i1,1) - x(i4,1);
    const F_FLOAT delz = x(i1,2) - x(i4,2);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const F_FLOAT r2inv = 1.0/rsq;
    const F_FLOAT r6inv = r2inv*r2inv*r2inv;

    F_FLOAT forcecoul;
    if (implicit) forcecoul = qqrd2e * q[i1]*q[i4]*r2inv;
    else forcecoul = qqrd2e * q[i1]*q[i4]*sqrt(r2inv);
    const F_FLOAT forcelj = r6inv * (d_lj14_1(itype,jtype)*r6inv - d_lj14_2(itype,jtype));
    const F_FLOAT fpair = d_weight[type] * (forcelj+forcecoul)*r2inv;

    F_FLOAT ecoul = 0.0;
    F_FLOAT evdwl = 0.0;
    if (eflag) {
      ecoul = d_weight[type] * forcecoul;
      evdwl = r6inv * (d_lj14_3(itype,jtype)*r6inv - d_lj14_4(itype,jtype));
      evdwl *= d_weight[type];
    }

    if (newton_bond || i1 < nlocal) {
      a_f(i1,0) += delx*fpair;
      a_f(i1,1) += dely*fpair;
      a_f(i1,2) += delz*fpair;
    }
    if (newton_bond || i4 < nlocal) {
      a_f(i4,0) -= delx*fpair;
      a_f(i4,1) -= dely*fpair;
      a_f(i4,2) -= delz*fpair;
    }

    if (EVFLAG) ev_tally(evm,i1,i4,evdwl,ecoul,fpair,delx,dely,delz);
  }
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void DihedralCharmmKokkos<DeviceType>::operator()(TagDihedralCharmmCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EVM_FLOAT evm;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralCharmmCompute<NEWTON_BOND,EVFLAG>(), n, evm);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCharmmKokkos<DeviceType>::allocate()
{
  DihedralCharmm::allocate();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCharmmKokkos<DeviceType>::coeff(int narg, char **arg)
{
  DihedralCharmm::coeff(narg, arg);

  int nd = atom->ndihedraltypes;
  typename AT::tdual_ffloat_1d k_k("DihedralCharmm::k",nd+1);
  typename AT::tdual_ffloat_1d k_multiplicity("DihedralCharmm::multiplicity",nd+1);
  typename AT::tdual_ffloat_1d k_shift("DihedralCharmm::shift",nd+1);
  typename AT::tdual_ffloat_1d k_cos_shift("DihedralCharmm::cos_shift",nd+1);
  typename AT::tdual_ffloat_1d k_sin_shift("DihedralCharmm::sin_shift",nd+1);
  typename AT::tdual_ffloat_1d k_weight("DihedralCharmm::weight",nd+1);

  d_k = k_k.template view<DeviceType>();
  d_multiplicity = k_multiplicity.template view<DeviceType>();
  d_shift = k_shift.template view<DeviceType>();
  d_cos_shift = k_cos_shift.template view<DeviceType>();
  d_sin_shift = k_sin_shift.template view<DeviceType>();
  d_weight = k_weight.template view<DeviceType>();

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_multiplicity.h_view[i] = multiplicity[i];
    k_shift.h_view[i] = shift[i];
    k_cos_shift.h_view[i] = cos_shift[i];
    k_sin_shift.h_view[i] = sin_shift[i];
    k_weight.h_view[i] = weight[i];
  }

  k_k.template modify<LMPHostType>();
  k_multiplicity.template modify<LMPHostType>();
  k_shift.template modify<LMPHostType>();
  k_cos_shift.template modify<LMPHostType>();
  k_sin_shift.template modify<LMPHostType>();
  k_weight.template modify<LMPHostType>();

  k_k.template sync<DeviceType>();
  k_multiplicity.template sync<DeviceType>();
  k_shift.template sync<DeviceType>();
  k_cos_shift.template sync<DeviceType>();
  k_sin_shift.template sync<DeviceType>();
  k_weight.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   error check and initialize all values needed for force computation
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCharmmKokkos<DeviceType>::init_style()
{
  DihedralCharmm::init_style();

  int n = atom->ntypes;
  DAT::tdual_ffloat_2d k_lj14_1("DihedralCharmm:lj14_1",n+1,n+1);
  DAT::tdual_ffloat_2d k_lj14_2("DihedralCharmm:lj14_2",n+1,n+1);
  DAT::tdual_ffloat_2d k_lj14_3("DihedralCharmm:lj14_3",n+1,n+1);
  DAT::tdual_ffloat_2d k_lj14_4("DihedralCharmm:lj14_4",n+1,n+1);

  d_lj14_1 = k_lj14_1.template view<DeviceType>();
  d_lj14_2 = k_lj14_2.template view<DeviceType>();
  d_lj14_3 = k_lj14_3.template view<DeviceType>();
  d_lj14_4 = k_lj14_4.template view<DeviceType>();


  if (weightflag) {
    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        k_lj14_1.h_view(i,j) = lj14_1[i][j];
        k_lj14_2.h_view(i,j) = lj14_2[i][j];
        k_lj14_3.h_view(i,j) = lj14_3[i][j];
        k_lj14_4.h_view(i,j) = lj14_4[i][j];
      }
    }
  }

  k_lj14_1.template modify<LMPHostType>();
  k_lj14_2.template modify<LMPHostType>();
  k_lj14_3.template modify<LMPHostType>();
  k_lj14_4.template modify<LMPHostType>();

  k_lj14_1.template sync<DeviceType>();
  k_lj14_2.template sync<DeviceType>();
  k_lj14_3.template sync<DeviceType>();
  k_lj14_4.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCharmmKokkos<DeviceType>::read_restart(FILE *fp)
{
  DihedralCharmm::read_restart(fp);

  int nd = atom->ndihedraltypes;
  typename AT::tdual_ffloat_1d k_k("DihedralCharmm::k",nd+1);
  typename AT::tdual_ffloat_1d k_multiplicity("DihedralCharmm::multiplicity",nd+1);
  typename AT::tdual_ffloat_1d k_shift("DihedralCharmm::shift",nd+1);
  typename AT::tdual_ffloat_1d k_cos_shift("DihedralCharmm::cos_shift",nd+1);
  typename AT::tdual_ffloat_1d k_sin_shift("DihedralCharmm::sin_shift",nd+1);
  typename AT::tdual_ffloat_1d k_weight("DihedralCharmm::weight",nd+1);

  d_k = k_k.template view<DeviceType>();
  d_multiplicity = k_multiplicity.template view<DeviceType>();
  d_shift = k_shift.template view<DeviceType>();
  d_cos_shift = k_cos_shift.template view<DeviceType>();
  d_sin_shift = k_sin_shift.template view<DeviceType>();
  d_weight = k_weight.template view<DeviceType>();

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_multiplicity.h_view[i] = multiplicity[i];
    k_shift.h_view[i] = shift[i];
    k_cos_shift.h_view[i] = cos_shift[i];
    k_sin_shift.h_view[i] = sin_shift[i];
    k_weight.h_view[i] = weight[i];
  }

  k_k.template modify<LMPHostType>();
  k_multiplicity.template modify<LMPHostType>();
  k_shift.template modify<LMPHostType>();
  k_cos_shift.template modify<LMPHostType>();
  k_sin_shift.template modify<LMPHostType>();
  k_weight.template modify<LMPHostType>();

  k_k.template sync<DeviceType>();
  k_multiplicity.template sync<DeviceType>();
  k_shift.template sync<DeviceType>();
  k_cos_shift.template sync<DeviceType>();
  k_sin_shift.template sync<DeviceType>();
  k_weight.template sync<DeviceType>();
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
void DihedralCharmmKokkos<DeviceType>::ev_tally(EVM_FLOAT &evm, const int i1, const int i2, const int i3, const int i4,
                        F_FLOAT &edihedral, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                        const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                        const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                        const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const
{
  E_FLOAT edihedralquarter;
  F_FLOAT v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) evm.emol += edihedral;
      else {
        edihedralquarter = 0.25*edihedral;
        if (i1 < nlocal) evm.emol += edihedralquarter;
        if (i2 < nlocal) evm.emol += edihedralquarter;
        if (i3 < nlocal) evm.emol += edihedralquarter;
        if (i4 < nlocal) evm.emol += edihedralquarter;
      }
    }
    if (eflag_atom) {
      edihedralquarter = 0.25*edihedral;
      if (newton_bond || i1 < nlocal) d_eatom[i1] += edihedralquarter;
      if (newton_bond || i2 < nlocal) d_eatom[i2] += edihedralquarter;
      if (newton_bond || i3 < nlocal) d_eatom[i3] += edihedralquarter;
      if (newton_bond || i4 < nlocal) d_eatom[i4] += edihedralquarter;
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
        evm.v[0] += v[0];
        evm.v[1] += v[1];
        evm.v[2] += v[2];
        evm.v[3] += v[3];
        evm.v[4] += v[4];
        evm.v[5] += v[5];
      } else {
        if (i1 < nlocal) {
          evm.v[0] += 0.25*v[0];
          evm.v[1] += 0.25*v[1];
          evm.v[2] += 0.25*v[2];
          evm.v[3] += 0.25*v[3];
          evm.v[4] += 0.25*v[4];
          evm.v[5] += 0.25*v[5];
        }
        if (i2 < nlocal) {
          evm.v[0] += 0.25*v[0];
          evm.v[1] += 0.25*v[1];
          evm.v[2] += 0.25*v[2];
          evm.v[3] += 0.25*v[3];
          evm.v[4] += 0.25*v[4];
          evm.v[5] += 0.25*v[5];
        }
        if (i3 < nlocal) {
          evm.v[0] += 0.25*v[0];
          evm.v[1] += 0.25*v[1];
          evm.v[2] += 0.25*v[2];
          evm.v[3] += 0.25*v[3];
          evm.v[4] += 0.25*v[4];
          evm.v[5] += 0.25*v[5];
        }
        if (i4 < nlocal) {
          evm.v[0] += 0.25*v[0];
          evm.v[1] += 0.25*v[1];
          evm.v[2] += 0.25*v[2];
          evm.v[3] += 0.25*v[3];
          evm.v[4] += 0.25*v[4];
          evm.v[5] += 0.25*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i1 < nlocal) {
        d_vatom(i1,0) += 0.25*v[0];
        d_vatom(i1,1) += 0.25*v[1];
        d_vatom(i1,2) += 0.25*v[2];
        d_vatom(i1,3) += 0.25*v[3];
        d_vatom(i1,4) += 0.25*v[4];
        d_vatom(i1,5) += 0.25*v[5];
      }
      if (newton_bond || i2 < nlocal) {
        d_vatom(i2,0) += 0.25*v[0];
        d_vatom(i2,1) += 0.25*v[1];
        d_vatom(i2,2) += 0.25*v[2];
        d_vatom(i2,3) += 0.25*v[3];
        d_vatom(i2,4) += 0.25*v[4];
        d_vatom(i2,5) += 0.25*v[5];
      }
      if (newton_bond || i3 < nlocal) {
        d_vatom(i3,0) += 0.25*v[0];
        d_vatom(i3,1) += 0.25*v[1];
        d_vatom(i3,2) += 0.25*v[2];
        d_vatom(i3,3) += 0.25*v[3];
        d_vatom(i3,4) += 0.25*v[4];
        d_vatom(i3,5) += 0.25*v[5];
      }
      if (newton_bond || i4 < nlocal) {
        d_vatom(i4,0) += 0.25*v[0];
        d_vatom(i4,1) += 0.25*v[1];
        d_vatom(i4,2) += 0.25*v[2];
        d_vatom(i4,3) += 0.25*v[3];
        d_vatom(i4,4) += 0.25*v[4];
        d_vatom(i4,5) += 0.25*v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void DihedralCharmmKokkos<DeviceType>::ev_tally(EVM_FLOAT &evm, const int i, const int j,
      const F_FLOAT &evdwl, const F_FLOAT &ecoul, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  E_FLOAT evdwlhalf,ecoulhalf,epairhalf;
  F_FLOAT v[6];


  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) {
        evm.evdwl += evdwl;
        evm.ecoul += ecoul;
      } else {
        evdwlhalf = 0.5*evdwl;
        ecoulhalf = 0.5*ecoul;
        if (i < nlocal) {
          evm.evdwl += evdwlhalf;
          evm.ecoul += ecoulhalf;
        }
        if (j < nlocal) {
          evm.evdwl += evdwlhalf;
          evm.ecoul += ecoulhalf;
        }
      }
    }
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_bond || i < nlocal) d_eatom_pair[i] += epairhalf;
      if (newton_bond || j < nlocal) d_eatom_pair[j] += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (vflag_global) {
      if (newton_bond) {
        evm.vp[0] += v[0];
        evm.vp[1] += v[1];
        evm.vp[2] += v[2];
        evm.vp[3] += v[3];
        evm.vp[4] += v[4];
        evm.vp[5] += v[5];
      } else {
        if (i < nlocal) {
          evm.vp[0] += 0.5*v[0];
          evm.vp[1] += 0.5*v[1];
          evm.vp[2] += 0.5*v[2];
          evm.vp[3] += 0.5*v[3];
          evm.vp[4] += 0.5*v[4];
          evm.vp[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          evm.vp[0] += 0.5*v[0];
          evm.vp[1] += 0.5*v[1];
          evm.vp[2] += 0.5*v[2];
          evm.vp[3] += 0.5*v[3];
          evm.vp[4] += 0.5*v[4];
          evm.vp[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        d_vatom_pair(i,0) += 0.5*v[0];
        d_vatom_pair(i,1) += 0.5*v[1];
        d_vatom_pair(i,2) += 0.5*v[2];
        d_vatom_pair(i,3) += 0.5*v[3];
        d_vatom_pair(i,4) += 0.5*v[4];
        d_vatom_pair(i,5) += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        d_vatom_pair(j,0) += 0.5*v[0];
        d_vatom_pair(j,1) += 0.5*v[1];
        d_vatom_pair(j,2) += 0.5*v[2];
        d_vatom_pair(j,3) += 0.5*v[3];
        d_vatom_pair(j,4) += 0.5*v[4];
        d_vatom_pair(j,5) += 0.5*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class DihedralCharmmKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class DihedralCharmmKokkos<LMPHostType>;
#endif
}

