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
   Contributing author: Mitch Murphy (alphataubio@gmail.com)
------------------------------------------------------------------------- */

#include "angle_spica_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "respa.h"
#include "update.h"

#include "lj_spica_common.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace LJSPICAParms;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleSPICAKokkos<DeviceType>::AngleSPICAKokkos(LAMMPS *lmp) : AngleSPICA(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleSPICAKokkos<DeviceType>::~AngleSPICAKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleSPICAKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"angle:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"angle:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  k_k.template sync<DeviceType>();
  k_theta0.template sync<DeviceType>();
  k_repscale.template sync<DeviceType>();
  k_lj_type.template sync<DeviceType>();
  k_lj1.template sync<DeviceType>();
  k_lj2.template sync<DeviceType>();
  k_lj3.template sync<DeviceType>();
  k_lj4.template sync<DeviceType>();
  k_rminsq.template sync<DeviceType>();
  k_emin.template sync<DeviceType>();


  // "It has to do with overlapping host/device in verlet_kokkos.cpp. For this reason, all topology styles (bond, angle, etc.) must set DATAMASK_READ, DATAMASK_MODIFY in the constructor and must not use atomKK->sync/modified. This is a gotcha that needed to be better documented."
  // https://matsci.org/t/a-few-kokkos-development-questions/56598
  //
  // atomKK->sync(execution_space,datamask_read);
  // if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  // else atomKK->modified(execution_space,F_MASK);
  //atomKK->k_type.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_anglelist.template sync<DeviceType>();
  anglelist = neighborKK->k_anglelist.template view<DeviceType>();
  int nanglelist = neighborKK->nanglelist;
  d_type = atomKK->k_type.template view<DeviceType>();
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleSPICACompute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleSPICACompute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleSPICACompute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleSPICACompute<0,0> >(0,nanglelist),*this);
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

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleSPICAKokkos<DeviceType>::operator()(TagAngleSPICACompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = anglelist(n,0);
  const int i2 = anglelist(n,1);
  const int i3 = anglelist(n,2);
  const int type = anglelist(n,3);

  // 1st bond

  const F_FLOAT delx1 = x(i1,0) - x(i2,0);
  const F_FLOAT dely1 = x(i1,1) - x(i2,1);
  const F_FLOAT delz1 = x(i1,2) - x(i2,2);

  const F_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  const F_FLOAT r1 = sqrt(rsq1);

  // 2nd bond

  const F_FLOAT delx2 = x(i3,0) - x(i2,0);
  const F_FLOAT dely2 = x(i3,1) - x(i2,1);
  const F_FLOAT delz2 = x(i3,2) - x(i2,2);

  const F_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  const F_FLOAT r2 = sqrt(rsq2);

  // angle (cos and sin)

  F_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  F_FLOAT s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  // 1-3 LJ interaction.
  // we only want to use the repulsive part,
  // and it can be scaled (or off).
  // so this has to be done here and not in the
  // general non-bonded code.

  F_FLOAT f13, e13, delx3, dely3, delz3;
  f13 = e13 = delx3 = dely3 = delz3 = 0.0;

  if (repflag) {

    delx3 = x(i1,0) - x(i3,0);
    dely3 = x(i1,1) - x(i3,1);
    delz3 = x(i1,2) - x(i3,2);
    const F_FLOAT rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;

    const int type1 = d_type[i1];
    const int type3 = d_type[i3];

    f13=0.0;
    e13=0.0;

    if (rsq3 < d_rminsq(type1,type3)) {
      const int ljt = d_lj_type(type1,type3);
      const double r2inv = 1.0/rsq3;

      if (ljt == LJ12_4) {
        const double r4inv=r2inv*r2inv;

        f13 = r4inv*(d_lj1(type1,type3)*r4inv*r4inv - d_lj2(type1,type3));
        if (eflag) e13 = r4inv*(d_lj3(type1,type3)*r4inv*r4inv - d_lj4(type1,type3));

      } else if (ljt == LJ9_6) {
        const double r3inv = r2inv*sqrt(r2inv);
        const double r6inv = r3inv*r3inv;

        f13 = r6inv*(d_lj1(type1,type3)*r3inv - d_lj2(type1,type3));
        if (eflag) e13 = r6inv*(d_lj3(type1,type3)*r3inv - d_lj4(type1,type3));

      } else if (ljt == LJ12_6) {
        const double r6inv = r2inv*r2inv*r2inv;

        f13 = r6inv*(d_lj1(type1,type3)*r6inv - d_lj2(type1,type3));
        if (eflag) e13 = r6inv*(d_lj3(type1,type3)*r6inv - d_lj4(type1,type3));

      } else if (ljt == LJ12_5) {
        const double r5inv = r2inv*r2inv*sqrt(r2inv);
        const double r7inv = r5inv*r2inv;

        f13 = r5inv*(d_lj1(type1,type3)*r7inv - d_lj2(type1,type3));
        if (eflag) e13 = r5inv*(d_lj3(type1,type3)*r7inv - d_lj4(type1,type3));
      }

      // make sure energy is 0.0 at the cutoff.
      if (eflag) e13 -= d_emin(type1,type3);

      f13 *= r2inv;
    }
  }

  // force & energy

  const F_FLOAT dtheta = acos(c) - d_theta0[type];
  const F_FLOAT tk = d_k[type] * dtheta;

  F_FLOAT eangle = 0.0;
  if (eflag) eangle = tk*dtheta;

  const F_FLOAT a = -2.0 * tk * s;
  const F_FLOAT a11 = a*c / rsq1;
  const F_FLOAT a12 = -a / (r1*r2);
  const F_FLOAT a22 = a*c / rsq2;

  F_FLOAT f1[3],f3[3];
  f1[0] = a11*delx1 + a12*delx2;
  f1[1] = a11*dely1 + a12*dely2;
  f1[2] = a11*delz1 + a12*delz2;
  f3[0] = a22*delx2 + a12*delx1;
  f3[1] = a22*dely2 + a12*dely1;
  f3[2] = a22*delz2 + a12*delz1;

  // apply force to each of 3 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += f1[0] + f13*delx3;
    a_f(i1,1) += f1[1] + f13*dely3;
    a_f(i1,2) += f1[2] + f13*delz3;
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) -= f1[0] + f3[0];
    a_f(i2,1) -= f1[1] + f3[1];
    a_f(i2,2) -= f1[2] + f3[2];
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += f3[0] - f13*delx3;
    a_f(i3,1) += f3[1] - f13*dely3;
    a_f(i3,2) += f3[2] - f13*delz3;
  }

  if (EVFLAG) {
    ev_tally(ev,i1,i2,i3,eangle,f1,f3,delx1,dely1,delz1,delx2,dely2,delz2);

    if (repflag)
      ev_tally13(ev,i1,i3,e13,f13,delx3,dely3,delz3);
  }
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleSPICAKokkos<DeviceType>::operator()(TagAngleSPICACompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleSPICACompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleSPICAKokkos<DeviceType>::allocate()
{
  AngleSPICA::allocate();

  int nangletypes = atom->nangletypes;
  k_k = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleSPICA::k",nangletypes+1);
  k_theta0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleSPICA::theta0",nangletypes+1);
  k_repscale = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleSPICA::repscale",nangletypes+1);
  k_setflag = typename ArrayTypes<DeviceType>::tdual_int_1d("AngleSPICA::setflag",nangletypes+1);

  d_k = k_k.template view<DeviceType>();
  d_theta0 = k_theta0.template view<DeviceType>();
  d_repscale = k_repscale.template view<DeviceType>();
  d_setflag = k_setflag.template view<DeviceType>();

  int ntypes = atom->ntypes;
  k_lj_type = typename ArrayTypes<DeviceType>::tdual_int_2d("AngleSPICA::lj_type",ntypes+1,ntypes+1);
  k_lj1 = typename ArrayTypes<DeviceType>::tdual_ffloat_2d("AngleSPICA::lj1",ntypes+1,ntypes+1);
  k_lj2 = typename ArrayTypes<DeviceType>::tdual_ffloat_2d("AngleSPICA::lj2",ntypes+1,ntypes+1);
  k_lj3 = typename ArrayTypes<DeviceType>::tdual_ffloat_2d("AngleSPICA::lj3",ntypes+1,ntypes+1);
  k_lj4 = typename ArrayTypes<DeviceType>::tdual_ffloat_2d("AngleSPICA::lj4",ntypes+1,ntypes+1);
  k_rminsq = typename ArrayTypes<DeviceType>::tdual_ffloat_2d("AngleSPICA::rminsq",ntypes+1,ntypes+1);
  k_emin = typename ArrayTypes<DeviceType>::tdual_ffloat_2d("AngleSPICA::emin",ntypes+1,ntypes+1);

  d_lj_type = k_lj_type.template view<DeviceType>();
  d_lj1 = k_lj1.template view<DeviceType>();
  d_lj2 = k_lj2.template view<DeviceType>();
  d_lj3 = k_lj3.template view<DeviceType>();
  d_lj4 = k_lj4.template view<DeviceType>();
  d_rminsq = k_rminsq.template view<DeviceType>();
  d_emin = k_emin.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleSPICAKokkos<DeviceType>::init_style()
{
  AngleSPICA::init_style();

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++) {
    for (int j = 1; j <= ntypes; j++) {
      k_lj_type.h_view(i,j) = lj_type[i][j];
      k_lj1.h_view(i,j) = lj1[i][j];
      k_lj2.h_view(i,j) = lj2[i][j];
      k_lj3.h_view(i,j) = lj3[i][j];
      k_lj4.h_view(i,j) = lj4[i][j];
      k_rminsq.h_view(i,j) = rminsq[i][j];
      k_emin.h_view(i,j) = emin[i][j];
    }
  }

  k_lj_type.template modify<LMPHostType>();
  k_lj1.template modify<LMPHostType>();
  k_lj2.template modify<LMPHostType>();
  k_lj3.template modify<LMPHostType>();
  k_lj4.template modify<LMPHostType>();
  k_rminsq.template modify<LMPHostType>();
  k_emin.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleSPICAKokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleSPICA::coeff(narg, arg);

  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_theta0.h_view[i] = theta0[i];
    k_repscale.h_view[i] = repscale[i];
    k_setflag.h_view[i] = setflag[i];
  }

  k_k.template modify<LMPHostType>();
  k_theta0.template modify<LMPHostType>();
  k_repscale.template modify<LMPHostType>();
  k_setflag.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleSPICAKokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleSPICA::read_restart(fp);

  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_theta0.h_view[i] = theta0[i];
    k_repscale.h_view[i] = repscale[i];
    k_setflag.h_view[i] = setflag[i];
  }

  k_k.template modify<LMPHostType>();
  k_theta0.template modify<LMPHostType>();
  k_repscale.template modify<LMPHostType>();
  k_setflag.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
KOKKOS_INLINE_FUNCTION
void AngleSPICAKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     F_FLOAT &eangle, F_FLOAT *f1, F_FLOAT *f3,
                     const F_FLOAT &delx1, const F_FLOAT &dely1, const F_FLOAT &delz1,
                     const F_FLOAT &delx2, const F_FLOAT &dely2, const F_FLOAT &delz2) const
{
  E_FLOAT eanglethird;
  F_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.template view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.template view<DeviceType>();

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += eangle;
      else {
        eanglethird = THIRD*eangle;

        if (i < nlocal) ev.evdwl += eanglethird;
        if (j < nlocal) ev.evdwl += eanglethird;
        if (k < nlocal) ev.evdwl += eanglethird;
      }
    }
    if (eflag_atom) {
      eanglethird = THIRD*eangle;

      if (newton_bond || i < nlocal) v_eatom[i] += eanglethird;
      if (newton_bond || j < nlocal) v_eatom[j] += eanglethird;
      if (newton_bond || k < nlocal) v_eatom[k] += eanglethird;
    }
  }

  if (vflag_either) {
    v[0] = delx1*f1[0] + delx2*f3[0];
    v[1] = dely1*f1[1] + dely2*f3[1];
    v[2] = delz1*f1[2] + delz2*f3[2];
    v[3] = delx1*f1[1] + delx2*f3[1];
    v[4] = delx1*f1[2] + delx2*f3[2];
    v[5] = dely1*f1[2] + dely2*f3[2];

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i < nlocal) {
          ev.v[0] += THIRD*v[0];
          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
        if (j < nlocal) {
          ev.v[0] += THIRD*v[0];
          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
        if (k < nlocal) {
          ev.v[0] += THIRD*v[0];

          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        v_vatom(i,0) += THIRD*v[0];
        v_vatom(i,1) += THIRD*v[1];
        v_vatom(i,2) += THIRD*v[2];
        v_vatom(i,3) += THIRD*v[3];
        v_vatom(i,4) += THIRD*v[4];
        v_vatom(i,5) += THIRD*v[5];
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j,0) += THIRD*v[0];
        v_vatom(j,1) += THIRD*v[1];
        v_vatom(j,2) += THIRD*v[2];
        v_vatom(j,3) += THIRD*v[3];
        v_vatom(j,4) += THIRD*v[4];
        v_vatom(j,5) += THIRD*v[5];
      }
      if (newton_bond || k < nlocal) {
        v_vatom(k,0) += THIRD*v[0];
        v_vatom(k,1) += THIRD*v[1];
        v_vatom(k,2) += THIRD*v[2];
        v_vatom(k,3) += THIRD*v[3];
        v_vatom(k,4) += THIRD*v[4];
        v_vatom(k,5) += THIRD*v[5];

      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void AngleSPICAKokkos<DeviceType>::ev_tally13(EV_FLOAT &ev, const int i, const int j,
                     const F_FLOAT &evdwl, const F_FLOAT &fpair,
                     const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const
{
  double v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.template view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.template view<DeviceType>();

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) {
        ev.evdwl += evdwl;
      } else {
        if (i < nlocal)
          ev.evdwl += 0.5*evdwl;
        if (j < nlocal)
          ev.evdwl += 0.5*evdwl;
      }
    }
    if (eflag_atom) {
      if (newton_bond || i < nlocal) v_eatom[i] += 0.5*evdwl;
      if (newton_bond || j < nlocal) v_eatom[j] += 0.5*evdwl;
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
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i < nlocal) {
          ev.v[0] += 0.5*v[0];
          ev.v[1] += 0.5*v[1];
          ev.v[2] += 0.5*v[2];
          ev.v[3] += 0.5*v[3];
          ev.v[4] += 0.5*v[4];
          ev.v[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          ev.v[0] += 0.5*v[0];
          ev.v[1] += 0.5*v[1];
          ev.v[2] += 0.5*v[2];
          ev.v[3] += 0.5*v[3];
          ev.v[4] += 0.5*v[4];
          ev.v[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        v_vatom(i,0) += 0.5*v[0];
        v_vatom(i,1) += 0.5*v[1];
        v_vatom(i,2) += 0.5*v[2];
        v_vatom(i,3) += 0.5*v[3];
        v_vatom(i,4) += 0.5*v[4];
        v_vatom(i,5) += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j,0) += 0.5*v[0];
        v_vatom(j,1) += 0.5*v[1];
        v_vatom(j,2) += 0.5*v[2];
        v_vatom(j,3) += 0.5*v[3];
        v_vatom(j,4) += 0.5*v[4];
        v_vatom(j,5) += 0.5*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class AngleSPICAKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleSPICAKokkos<LMPHostType>;
#endif
}

