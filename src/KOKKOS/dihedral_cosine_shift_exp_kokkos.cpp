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
   [ based on dihedral_harmonic_kokkos.cpp and dihedral_cosine_shift_exp.cpp ]
------------------------------------------------------------------------- */

#include "dihedral_cosine_shift_exp_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr double TOLERANCE = 0.05;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralCosineShiftExpKokkos<DeviceType>::DihedralCosineShiftExpKokkos(LAMMPS *lmp)
  : DihedralCosineShiftExp(lmp)
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
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralCosineShiftExpKokkos<DeviceType>::~DihedralCosineShiftExpKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCosineShiftExpKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  k_umin.template sync<DeviceType>();
  k_a.template sync<DeviceType>();
  k_opt1.template sync<DeviceType>();
  k_cost.template sync<DeviceType>();
  k_sint.template sync<DeviceType>();
  k_doExpansion.template sync<DeviceType>();

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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralCosineShiftExpCompute<1,1> >(0,ndihedrallist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralCosineShiftExpCompute<0,1> >(0,ndihedrallist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralCosineShiftExpCompute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralCosineShiftExpCompute<0,0> >(0,ndihedrallist),*this);
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

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralCosineShiftExpKokkos<DeviceType>::operator()(TagDihedralCosineShiftExpCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

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

  // c,s calculation

  const KK_FLOAT ax = vb1y*vb2zm - vb1z*vb2ym;
  const KK_FLOAT ay = vb1z*vb2xm - vb1x*vb2zm;
  const KK_FLOAT az = vb1x*vb2ym - vb1y*vb2xm;
  const KK_FLOAT bx = vb3y*vb2zm - vb3z*vb2ym;
  const KK_FLOAT by = vb3z*vb2xm - vb3x*vb2zm;
  const KK_FLOAT bz = vb3x*vb2ym - vb3y*vb2xm;

  const KK_FLOAT rasq = ax*ax + ay*ay + az*az;
  const KK_FLOAT rbsq = bx*bx + by*by + bz*bz;
  const KK_FLOAT rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
  const KK_FLOAT rg = sqrt(rgsq);

  KK_FLOAT rginv,ra2inv,rb2inv;
  rginv = ra2inv = rb2inv = 0.0;
  if (rg > 0) rginv = 1.0/rg;
  if (rasq > 0) ra2inv = 1.0/rasq;
  if (rbsq > 0) rb2inv = 1.0/rbsq;
  const KK_FLOAT rabinv = sqrt(ra2inv*rb2inv);

  KK_FLOAT c = (ax*bx + ay*by + az*bz)*rabinv;
  const KK_FLOAT s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

  // error check

  if ((c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) && !d_warning_flag())
    d_warning_flag() = 1;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  const KK_FLOAT aa = d_a[type];
  const KK_FLOAT uumin = d_umin[type];

  const KK_FLOAT cccpsss = c*d_cost[type] + s*d_sint[type];
  const KK_FLOAT cssmscc = c*d_sint[type] - s*d_cost[type];

  KK_FLOAT edihedral = 0.0;
  KK_FLOAT df;

  if (d_doExpansion[type]) {
    // |a|<0.001 so use expansions, relative precision <1e-5
    if (EVFLAG && eflag) edihedral = -0.125*(1+cccpsss)*(4+aa*(cccpsss-1))*uumin;
    df = 0.5*uumin*(cssmscc + 0.5*aa*cccpsss);
  } else {
    const KK_FLOAT exp2 = exp(0.5*aa*(1+cccpsss));
    if (EVFLAG && eflag) edihedral = d_opt1[type]*(1-exp2);
    df = 0.5*d_opt1[type]*aa*(exp2*cssmscc);
  }

  const KK_FLOAT fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
  const KK_FLOAT hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
  const KK_FLOAT fga = fg*ra2inv*rginv;
  const KK_FLOAT hgb = hg*rb2inv*rginv;
  const KK_FLOAT gaa = -ra2inv*rg;
  const KK_FLOAT gbb = rb2inv*rg;

  const KK_FLOAT dtfx = gaa*ax;
  const KK_FLOAT dtfy = gaa*ay;
  const KK_FLOAT dtfz = gaa*az;
  const KK_FLOAT dtgx = fga*ax - hgb*bx;
  const KK_FLOAT dtgy = fga*ay - hgb*by;
  const KK_FLOAT dtgz = fga*az - hgb*bz;
  const KK_FLOAT dthx = gbb*bx;
  const KK_FLOAT dthy = gbb*by;
  const KK_FLOAT dthz = gbb*bz;

  const KK_FLOAT sx2 = df*dtgx;
  const KK_FLOAT sy2 = df*dtgy;
  const KK_FLOAT sz2 = df*dtgz;

  KK_FLOAT f1[3],f2[3],f3[3],f4[3];
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
    ev_tally(ev,i1,i2,i3,i4,edihedral,f1,f3,f4,
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void DihedralCosineShiftExpKokkos<DeviceType>::operator()(TagDihedralCosineShiftExpCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralCosineShiftExpCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCosineShiftExpKokkos<DeviceType>::allocate()
{
  DihedralCosineShiftExp::allocate();

  int n = atom->ndihedraltypes;
  k_umin = DAT::tdual_kkfloat_1d("DihedralCosineShiftExp::umin",n+1);
  k_a = DAT::tdual_kkfloat_1d("DihedralCosineShiftExp::a",n+1);
  k_opt1 = DAT::tdual_kkfloat_1d("DihedralCosineShiftExp::opt1",n+1);
  k_cost = DAT::tdual_kkfloat_1d("DihedralCosineShiftExp::cost",n+1);
  k_sint = DAT::tdual_kkfloat_1d("DihedralCosineShiftExp::sint",n+1);
  k_doExpansion = DAT::tdual_int_1d("DihedralCosineShiftExp::doExpansion",n+1);

  d_umin = k_umin.template view<DeviceType>();
  d_a = k_a.template view<DeviceType>();
  d_opt1 = k_opt1.template view<DeviceType>();
  d_cost = k_cost.template view<DeviceType>();
  d_sint = k_sint.template view<DeviceType>();
  d_doExpansion = k_doExpansion.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCosineShiftExpKokkos<DeviceType>::coeff(int narg, char **arg)
{
  DihedralCosineShiftExp::coeff(narg, arg);

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->ndihedraltypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_umin.view_host()[i] = umin[i];
    k_a.view_host()[i] = a[i];
    k_opt1.view_host()[i] = opt1[i];
    k_cost.view_host()[i] = cost[i];
    k_sint.view_host()[i] = sint[i];
    k_doExpansion.view_host()[i] = (int)doExpansion[i];
  }

  k_umin.modify_host();
  k_a.modify_host();
  k_opt1.modify_host();
  k_cost.modify_host();
  k_sint.modify_host();
  k_doExpansion.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralCosineShiftExpKokkos<DeviceType>::read_restart(FILE *fp)
{
  DihedralCosineShiftExp::read_restart(fp);

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_umin.view_host()[i] = umin[i];
    k_a.view_host()[i] = a[i];
    k_opt1.view_host()[i] = opt1[i];
    k_cost.view_host()[i] = cost[i];
    k_sint.view_host()[i] = sint[i];
    k_doExpansion.view_host()[i] = (int)doExpansion[i];
  }

  k_umin.modify_host();
  k_a.modify_host();
  k_opt1.modify_host();
  k_cost.modify_host();
  k_sint.modify_host();
  k_doExpansion.modify_host();
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
void DihedralCosineShiftExpKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
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
template class DihedralCosineShiftExpKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class DihedralCosineShiftExpKokkos<LMPHostType>;
#endif
}
