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

#include "pair_brownian_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "fix.h"
#include "fix_wall.h"
#include "force.h"
#include "input.h"
#include "kokkos.h"
#include "math_const.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecialKokkos;

// same as fix_wall.cpp

enum { EDGE, CONSTANT, VARIABLE };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBrownianKokkos<DeviceType>::PairBrownianKokkos(LAMMPS *lmp) : PairBrownian(lmp),
                                                                  rand_pool(seed + comm->me)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | VIRIAL_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBrownianKokkos<DeviceType>::~PairBrownianKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cut_inner,cut_inner);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBrownianKokkos<DeviceType>::init_style()
{
  PairBrownian::init_style();

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBrownianKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF)                             // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2) {    // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (int j = 0; j < 3; j++) dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)) {
        double wallhi[3], walllo[3];
        for (int j = 0; j < 3; j++) {
          wallhi[j] = domain->prd[j];
          walllo[j] = 0;
        }
        for (int m = 0; m < wallfix->nwall; m++) {
          int dim = wallfix->wallwhich[m] / 2;
          int side = wallfix->wallwhich[m] % 2;
          if (wallfix->xstyle[m] == VARIABLE) {
            wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
          } else
            wallcoord = wallfix->coord0[m];
          if (side == 0)
            walllo[dim] = wallcoord;
          else
            wallhi[dim] = wallcoord;
        }
        for (int j = 0; j < 3; j++) dims[j] = wallhi[j] - walllo[j];
      }
      double vol_T = dims[0] * dims[1] * dims[2];
      double vol_f = vol_P / vol_T;
      if (flaglog == 0) {
        R0 = 6 * MY_PI * mu * rad * (1.0 + 2.16 * vol_f);
        RT0 = 8 * MY_PI * mu * cube(rad);
        //RS0 = 20.0/3.0*MY_PI*mu*pow(rad,3)*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
      } else {
        R0 = 6 * MY_PI * mu * rad * (1.0 + 2.725 * vol_f - 6.583 * vol_f * vol_f);
        RT0 = 8 * MY_PI * mu * cube(rad) * (1.0 + 0.749 * vol_f - 2.469 * vol_f * vol_f);
        //RS0 = 20.0/3.0*MY_PI*mu*pow(rad,3)*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
      }
    }

  // scale factor for Brownian moments

  prethermostat = sqrt(24.0 * force->boltz * t_target / update->dt);
  prethermostat *= sqrt(force->vxmu2f / force->ftm2v / force->mvv2e);

  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_cut_inner.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  vxmu2f = force->vxmu2f;

  // loop over neighbors of my atoms

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  copymode = 1;

  EV_FLOAT ev;

  if (flagfld) { // FLAGFLD == 1
    if (vflag_either) { // VFLAG == 1
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,1,1,1> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,0,1,1> >(0,inum),*this,ev);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,1,1,1> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,0,1,1> >(0,inum),*this,ev);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,1,1,1> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,0,1,1> >(0,inum),*this,ev);
      }
    } else {  // VFLAG==0
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,1,0,1> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,0,0,1> >(0,inum),*this);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,1,0,1> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,0,0,1> >(0,inum),*this);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,1,0,1> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,0,0,1> >(0,inum),*this);
      }
    }
  } else { // FLAGFLD == 0
    if (evflag) { // VFLAG== 1
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,1,1,0> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,0,1,0> >(0,inum),*this,ev);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,1,1,0> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,0,1,0> >(0,inum),*this,ev);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,1,1,0> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,0,1,0> >(0,inum),*this,ev);
      }
    } else { // VFLAG == 0
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,1,0,0> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALF,0,0,0> >(0,inum),*this);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,1,0,0> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<HALFTHREAD,0,0,0> >(0,inum),*this);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,1,0,0> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBrownianCompute<FULL,0,0,0> >(0,inum),*this);
      }
    }
  }

  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
KOKKOS_INLINE_FUNCTION
void PairBrownianKokkos<DeviceType>::operator()(TagPairBrownianCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  rand_type rand_gen = rand_pool.get_state();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type[i];
  const LMP_FLOAT radi = radius[i];
  const int jnum = d_numneigh[i];

  LMP_FLOAT a_sq, a_sh, a_pu;
  LMP_FLOAT xl[3], p1[3], p2[3], p3[3];

  F_FLOAT fx_i = 0.0;
  F_FLOAT fy_i = 0.0;
  F_FLOAT fz_i = 0.0;

  F_FLOAT torquex_i = 0.0;
  F_FLOAT torquey_i = 0.0;
  F_FLOAT torquez_i = 0.0;

  if (FLAGFLD) {
    fx_i = prethermostat * sqrt(R0) * (rand_gen.drand() - 0.5);
    fy_i = prethermostat * sqrt(R0) * (rand_gen.drand() - 0.5);
    fz_i = prethermostat * sqrt(R0) * (rand_gen.drand() - 0.5);
    if (flaglog) {
      torquex_i = prethermostat * sqrt(RT0) * (rand_gen.drand() - 0.5);
      torquey_i = prethermostat * sqrt(RT0) * (rand_gen.drand() - 0.5);
      torquez_i = prethermostat * sqrt(RT0) * (rand_gen.drand() - 0.5);
    }
  }

  if (flagHI) {

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      const X_FLOAT delx = xtmp - x(j,0);
      const X_FLOAT dely = ytmp - x(j,1);
      const X_FLOAT delz = ztmp - x(j,2);
      const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      const int jtype = type[j];

      if(rsq < d_cutsq(itype,jtype)) {

        const LMP_FLOAT r = sqrt(rsq);

        // scalar resistances a_sq and a_sh

        LMP_FLOAT h_sep = r - 2.0 * radi;

        // if less than minimum gap, use minimum gap instead

        if (r < d_cut_inner(itype,jtype)) h_sep = d_cut_inner(itype,jtype) - 2.0 * radi;

        // scale h_sep by radi

        h_sep = h_sep / radi;

        if (flaglog) {
          a_sq = 6.0 * MY_PI * mu * radi * (1.0 / 4.0 / h_sep + 9.0 / 40.0 * log(1.0 / h_sep));
          a_sh = 6.0 * MY_PI * mu * radi * (1.0 / 6.0 * log(1.0 / h_sep));
          a_pu = 8.0 * MY_PI * mu * cube(radi) * (3.0 / 160.0 * log(1.0 / h_sep));
        } else
          a_sq = 6.0 * MY_PI * mu * radi * (1.0 / 4.0 / h_sep);

        // generate the Pairwise Brownian Force: a_sq

        LMP_FLOAT Fbmag = prethermostat * sqrt(a_sq);

        // generate a random number

        LMP_FLOAT randr = rand_gen.drand() - 0.5;

        // contribution due to Brownian motion

        F_FLOAT fx = Fbmag * randr * delx / r;
        F_FLOAT fy = Fbmag * randr * dely / r;
        F_FLOAT fz = Fbmag * randr * delz / r;

        // add terms due to a_sh

        if (flaglog) {

          // generate two orthogonal vectors to the line of centers

          p1[0] = delx / r;
          p1[1] = dely / r;
          p1[2] = delz / r;
          set_3_orthogonal_vectors(p1, p2, p3);

          // magnitude

          Fbmag = prethermostat * sqrt(a_sh);

          // force in each of the two directions

          randr = rand_gen.drand() - 0.5;

          fx += Fbmag * randr * p2[0];
          fy += Fbmag * randr * p2[1];
          fz += Fbmag * randr * p2[2];

          randr = rand_gen.drand() - 0.5;

          fx += Fbmag * randr * p3[0];
          fy += Fbmag * randr * p3[1];
          fz += Fbmag * randr * p3[2];
        }

        // scale forces to appropriate units

        fx = vxmu2f * fx;
        fy = vxmu2f * fy;
        fz = vxmu2f * fz;

        // sum to total force

        fx_i -= fx;
        fy_i -= fy;
        fz_i -= fz;

        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
          a_f(j,0) += fx;
          a_f(j,1) += fy;
          a_f(j,2) += fz;
        }

        // torque due to the Brownian Force

        if (flaglog) {

          // location of the point of closest approach on I from its center

          xl[0] = -delx / r * radi;
          xl[1] = -dely / r * radi;
          xl[2] = -delz / r * radi;

          // torque = xl_cross_F

          F_FLOAT tx = xl[1] * fz - xl[2] * fy;
          F_FLOAT ty = xl[2] * fx - xl[0] * fz;
          F_FLOAT tz = xl[0] * fy - xl[1] * fx;

          // torque is same on both particles

          torquex_i -= tx;
          torquey_i -= ty;
          torquez_i -= tz;

          if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
            a_torque(j,0) -= tx;
            a_torque(j,1) -= ty;
            a_torque(j,2) -= tz;
          }

          // torque due to a_pu

          Fbmag = prethermostat * sqrt(a_pu);

          // force in each direction

          randr = rand_gen.drand() - 0.5;

          tx = Fbmag * randr * p2[0];
          ty = Fbmag * randr * p2[1];
          tz = Fbmag * randr * p2[2];

          randr = rand_gen.drand() - 0.5;

          tx += Fbmag * randr * p3[0];
          ty += Fbmag * randr * p3[1];
          tz += Fbmag * randr * p3[2];

          // torque has opposite sign on two particles

          torquex_i -= tx;
          torquey_i -= ty;
          torquez_i -= tz;

          if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
            a_torque(j,0) += tx;
            a_torque(j,1) += ty;
            a_torque(j,2) += tz;
          }
        }

        if (VFLAG)
          ev_tally_xyz<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, -fx, -fy, -fz, delx, dely, delz);
      }
    }

  } // if(flagHI)

  rand_pool.free_state(rand_gen);

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torquex_i;
  a_torque(i,1) += torquey_i;
  a_torque(i,2) += torquez_i;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
KOKKOS_INLINE_FUNCTION
void PairBrownianKokkos<DeviceType>::operator()(TagPairBrownianCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>(TagPairBrownianCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>(), ii, ev);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairBrownianKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT & ev, int i, int j,
                                                        F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
                                                        X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const
{
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

  const F_FLOAT v0 = delx*fx;
  const F_FLOAT v1 = dely*fy;
  const F_FLOAT v2 = delz*fz;
  const F_FLOAT v3 = delx*fy;
  const F_FLOAT v4 = delx*fz;
  const F_FLOAT v5 = dely*fz;

  if (vflag_global) {
    if (NEIGHFLAG != FULL) {
      if (NEWTON_PAIR) { // neigh half, newton on
        ev.v[0] += v0;
        ev.v[1] += v1;
        ev.v[2] += v2;
        ev.v[3] += v3;
        ev.v[4] += v4;
        ev.v[5] += v5;
      } else { // neigh half, newton off
        if (i < nlocal) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
        if (j < nlocal) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
      }
    } else { //neigh full
      ev.v[0] += 0.5*v0;
      ev.v[1] += 0.5*v1;
      ev.v[2] += 0.5*v2;
      ev.v[3] += 0.5*v3;
      ev.v[4] += 0.5*v4;
      ev.v[5] += 0.5*v5;
    }
  }

  if (vflag_atom) {

    if (NEIGHFLAG == FULL || NEWTON_PAIR || i < nlocal) {
      v_vatom(i,0) += 0.5*v0;
      v_vatom(i,1) += 0.5*v1;
      v_vatom(i,2) += 0.5*v2;
      v_vatom(i,3) += 0.5*v3;
      v_vatom(i,4) += 0.5*v4;
      v_vatom(i,5) += 0.5*v5;
    }
    if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) {
      v_vatom(j,0) += 0.5*v0;
      v_vatom(j,1) += 0.5*v1;
      v_vatom(j,2) += 0.5*v2;
      v_vatom(j,3) += 0.5*v3;
      v_vatom(j,4) += 0.5*v4;
      v_vatom(j,5) += 0.5*v5;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBrownianKokkos<DeviceType>::allocate()
{
  PairBrownian::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_inner);
  memoryKK->create_kokkos(k_cut_inner,cut_inner,n+1,n+1,"pair:cut_inner");
  d_cut_inner = k_cut_inner.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBrownianKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg != 7 && narg != 9) error->all(FLERR, "Illegal pair_style command");

  PairBrownian::settings(narg,arg);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairBrownianKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBrownian::init_one(i,j);
  double cutinnerm = cut_inner[i][j];

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();

  k_cut_inner.h_view(i,j) = k_cut_inner.h_view(j,i) = cutinnerm;
  k_cut_inner.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBrownianKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairBrownian::coeff(narg,arg);
}

namespace LAMMPS_NS {
template class PairBrownianKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBrownianKokkos<LMPHostType>;
#endif
}

