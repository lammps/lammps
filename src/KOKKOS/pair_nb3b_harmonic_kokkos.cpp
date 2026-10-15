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
------------------------------------------------------------------------- */

#include "pair_nb3b_harmonic_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairNb3bHarmonicKokkos<DeviceType>::PairNb3bHarmonicKokkos(LAMMPS *lmp) : PairNb3bHarmonic(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairNb3bHarmonicKokkos<DeviceType>::~PairNb3bHarmonicKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairNb3bHarmonicKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  nall = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  EV_FLOAT ev;
  EV_FLOAT ev_all;

  int chunk_size = 0;
  if (lmp->kokkos->threads_per_atom_set)
    chunk_size = lmp->kokkos->threads_per_atom;

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if (((int) d_neighbors_short.extent(1) < max_neighs) ||
      ((int) d_neighbors_short.extent(0) < ignum)) {
    d_neighbors_short = typename AT::t_int_2d_dl("Nb3bHarmonic::neighbors_short",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short.extent(0) < ignum)
    d_numneigh_short = typename AT::t_int_1d("Nb3bHarmonic::numneighs_short",ignum*1.2);

  if (chunk_size)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairNb3bHarmonicComputeShortNeigh>(0,inum,Kokkos::ChunkSize(chunk_size)), *this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairNb3bHarmonicComputeShortNeigh>(0,inum), *this);

  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairNb3bHarmonicCompute<HALF,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairNb3bHarmonicCompute<HALF,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairNb3bHarmonicCompute<HALF,0> >(0,inum),*this);
    }
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairNb3bHarmonicCompute<HALFTHREAD,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairNb3bHarmonicCompute<HALFTHREAD,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairNb3bHarmonicCompute<HALFTHREAD,0> >(0,inum),*this);
    }
    ev_all += ev;
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += static_cast<double>(ev_all.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev_all.v[0]);
    virial[1] += static_cast<double>(ev_all.v[1]);
    virial[2] += static_cast<double>(ev_all.v[2]);
    virial[3] += static_cast<double>(ev_all.v[3]);
    virial[4] += static_cast<double>(ev_all.v[4]);
    virial[5] += static_cast<double>(ev_all.v[5]);
  }

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f            = {};
    dup_eatom        = {};
    dup_vatom        = {};
  }
}

/* ----------------------------------------------------------------------
   short neighbor list: neighbors within the I-J cutoff, which is also
   the I-K cutoff when the neighbor is used as atom K
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairNb3bHarmonicKokkos<DeviceType>::operator()(TagPairNb3bHarmonicComputeShortNeigh, const int& ii) const {
    const int i = d_ilist[ii];
    const int itype = d_map[type[i]];
    const KK_FLOAT xtmp = x(i,0);
    const KK_FLOAT ytmp = x(i,1);
    const KK_FLOAT ztmp = x(i,2);

    const int jnum = d_numneigh[i];
    int inside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const int jtype = d_map[type[j]];

      const KK_FLOAT delx = xtmp - x(j,0);
      const KK_FLOAT dely = ytmp - x(j,1);
      const KK_FLOAT delz = ztmp - x(j,2);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      const int ijparam = d_elem3param(itype,jtype,jtype);
      if (rsq <= static_cast<KK_FLOAT>(d_params[ijparam].cutsq)) {
        d_neighbors_short(ii,inside) = j;
        inside++;
      }
    }
    d_numneigh_short(ii) = inside;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairNb3bHarmonicKokkos<DeviceType>::operator()(TagPairNb3bHarmonicCompute<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  KK_FLOAT delr1[3],delr2[3];
  KK_ACC_FLOAT fj[3],fk[3];
  KK_FLOAT evdwl = 0.0;

  const int i = d_ilist[ii];
  const int itype = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  const int jnum = d_numneigh_short[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  const int jnumm1 = jnum - 1;

  for (int jj = 0; jj < jnumm1; jj++) {
    int j = d_neighbors_short(ii,jj);
    const int jtype = d_map[type[j]];
    const int ijparam = d_elem3param(itype,jtype,jtype);
    delr1[0] = x(j,0) - xtmp;
    delr1[1] = x(j,1) - ytmp;
    delr1[2] = x(j,2) - ztmp;
    const KK_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

    KK_ACC_FLOAT fxtmpj = 0.0;
    KK_ACC_FLOAT fytmpj = 0.0;
    KK_ACC_FLOAT fztmpj = 0.0;

    for (int kk = jj+1; kk < jnum; kk++) {
      int k = d_neighbors_short(ii,kk);
      const int ktype = d_map[type[k]];
      const int ikparam = d_elem3param(itype,ktype,ktype);
      const int ijkparam = d_elem3param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const KK_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

      threebody_kk(d_params[ijparam],d_params[ikparam],d_params[ijkparam],
                   rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

      fxtmpi -= fj[0] + fk[0];
      fytmpi -= fj[1] + fk[1];
      fztmpi -= fj[2] + fk[2];
      fxtmpj += fj[0];
      fytmpj += fj[1];
      fztmpj += fj[2];
      a_f(k,0) += fk[0];
      a_f(k,1) += fk[1];
      a_f(k,2) += fk[2];

      if (EVFLAG) {
        if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
        if (vflag_either || eflag_atom) this->template ev_tally3<NEIGHFLAG>(ev,i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
    }

    a_f(j,0) += fxtmpj;
    a_f(j,1) += fytmpj;
    a_f(j,2) += fztmpj;
  }

  a_f(i,0) += fxtmpi;
  a_f(i,1) += fytmpi;
  a_f(i,2) += fztmpi;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairNb3bHarmonicKokkos<DeviceType>::operator()(TagPairNb3bHarmonicCompute<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairNb3bHarmonicCompute<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairNb3bHarmonicKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairNb3bHarmonic::coeff(narg,arg);

  // sync map

  int n = atom->ntypes;

  DAT::tdual_int_1d k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.view_host();

  for (int i = 1; i <= n; i++)
    h_map[i] = map[i];

  k_map.modify_host();
  k_map.template sync<DeviceType>();

  d_map = k_map.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairNb3bHarmonicKokkos<DeviceType>::init_style()
{
  PairNb3bHarmonic::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair {}/kk",
               (variant == SCREENED) ? "nb3b/screened" : "nb3b/harmonic");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairNb3bHarmonicKokkos<DeviceType>::setup_params()
{
  PairNb3bHarmonic::setup_params();

  // sync elem3param and params

  DAT::tdual_int_3d k_elem3param = DAT::tdual_int_3d("pair:elem3param",nelements,nelements,nelements);
  HAT::t_int_3d h_elem3param = k_elem3param.view_host();

  tdual_param_1d k_params = tdual_param_1d("pair:params",nparams);
  t_host_param_1d h_params = k_params.view_host();

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++)
        h_elem3param(i,j,k) = elem3param[i][j][k];

  for (int m = 0; m < nparams; m++)
    h_params[m] = params[m];

  k_elem3param.modify_host();
  k_elem3param.template sync<DeviceType>();
  k_params.modify_host();
  k_params.template sync<DeviceType>();

  d_elem3param = k_elem3param.template view<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   harmonic angle term, optionally multiplied by the screening function
   exp(-rij/rhoij - rik/rhoik) of the nb3b/screened variant
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairNb3bHarmonicKokkos<DeviceType>::threebody_kk(const Param& paramij, const Param& paramik, const Param& paramijk,
                       const KK_FLOAT& rsq1, const KK_FLOAT& rsq2,
                       KK_FLOAT *delr1, KK_FLOAT *delr2,
                       KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, const int& eflag, KK_FLOAT& eng) const
{
  const KK_FLOAT r1 = Kokkos::sqrt(rsq1);
  const KK_FLOAT r2 = Kokkos::sqrt(rsq2);

  // angle (cos and sin)

  KK_FLOAT c = delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2];
  c /= r1*r2;

  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);

  KK_FLOAT s = Kokkos::sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  if (s < static_cast<KK_FLOAT>(SMALL)) s = static_cast<KK_FLOAT>(SMALL);
  s = static_cast<KK_FLOAT>(1.0)/s;

  // force & energy

  const KK_FLOAT dtheta = Kokkos::acos(c) - static_cast<KK_FLOAT>(paramijk.theta0);

  if (variant == SCREENED) {
    const KK_FLOAT rho1inv = static_cast<KK_FLOAT>(paramij.invrho);
    const KK_FLOAT rho2inv = static_cast<KK_FLOAT>(paramik.invrho);
    const KK_FLOAT scr = Kokkos::exp(-r1*rho1inv - r2*rho2inv);
    const KK_FLOAT tk = static_cast<KK_FLOAT>(paramijk.k_theta) * dtheta * scr;
    const KK_FLOAT t00 = tk*dtheta;

    if (eflag) eng = t00;

    const KK_FLOAT a = static_cast<KK_FLOAT>(-2.0) * tk * s;
    const KK_FLOAT a11 = a*c / rsq1;
    const KK_FLOAT a12 = -a / (r1*r2);
    const KK_FLOAT a22 = a*c / rsq2;
    const KK_FLOAT ratio1 = rho1inv / r1;
    const KK_FLOAT ratio2 = rho2inv / r2;

    fj[0] = static_cast<KK_ACC_FLOAT>(a11*delr1[0] + a12*delr2[0] + t00*ratio1*delr1[0]);
    fj[1] = static_cast<KK_ACC_FLOAT>(a11*delr1[1] + a12*delr2[1] + t00*ratio1*delr1[1]);
    fj[2] = static_cast<KK_ACC_FLOAT>(a11*delr1[2] + a12*delr2[2] + t00*ratio1*delr1[2]);
    fk[0] = static_cast<KK_ACC_FLOAT>(a22*delr2[0] + a12*delr1[0] + t00*ratio2*delr2[0]);
    fk[1] = static_cast<KK_ACC_FLOAT>(a22*delr2[1] + a12*delr1[1] + t00*ratio2*delr2[1]);
    fk[2] = static_cast<KK_ACC_FLOAT>(a22*delr2[2] + a12*delr1[2] + t00*ratio2*delr2[2]);
  } else {
    const KK_FLOAT tk = static_cast<KK_FLOAT>(paramijk.k_theta) * dtheta;

    if (eflag) eng = tk*dtheta;

    const KK_FLOAT a = static_cast<KK_FLOAT>(-2.0) * tk * s;
    const KK_FLOAT a11 = a*c / rsq1;
    const KK_FLOAT a12 = -a / (r1*r2);
    const KK_FLOAT a22 = a*c / rsq2;

    fj[0] = static_cast<KK_ACC_FLOAT>(a11*delr1[0] + a12*delr2[0]);
    fj[1] = static_cast<KK_ACC_FLOAT>(a11*delr1[1] + a12*delr2[1]);
    fj[2] = static_cast<KK_ACC_FLOAT>(a11*delr1[2] + a12*delr2[2]);
    fk[0] = static_cast<KK_ACC_FLOAT>(a22*delr2[0] + a12*delr1[0]);
    fk[1] = static_cast<KK_ACC_FLOAT>(a22*delr2[1] + a12*delr1[1]);
    fk[2] = static_cast<KK_ACC_FLOAT>(a22*delr2[2] + a12*delr1[2]);
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairNb3bHarmonicKokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
          const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                     KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const
{
  KK_FLOAT epairthird,v[6];

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    epairthird = static_cast<KK_FLOAT>(THIRD) * (evdwl + ecoul);
    a_eatom[i] += static_cast<KK_ACC_FLOAT>(epairthird);
    a_eatom[j] += static_cast<KK_ACC_FLOAT>(epairthird);
    a_eatom[k] += static_cast<KK_ACC_FLOAT>(epairthird);
  }

  if (vflag_either) {
    v[0] = drji[0]*static_cast<KK_FLOAT>(fj[0]) + drki[0]*static_cast<KK_FLOAT>(fk[0]);
    v[1] = drji[1]*static_cast<KK_FLOAT>(fj[1]) + drki[1]*static_cast<KK_FLOAT>(fk[1]);
    v[2] = drji[2]*static_cast<KK_FLOAT>(fj[2]) + drki[2]*static_cast<KK_FLOAT>(fk[2]);
    v[3] = drji[0]*static_cast<KK_FLOAT>(fj[1]) + drki[0]*static_cast<KK_FLOAT>(fk[1]);
    v[4] = drji[0]*static_cast<KK_FLOAT>(fj[2]) + drki[0]*static_cast<KK_FLOAT>(fk[2]);
    v[5] = drji[1]*static_cast<KK_FLOAT>(fj[2]) + drki[1]*static_cast<KK_FLOAT>(fk[2]);

    if (vflag_global) {
      ev.v[0] += static_cast<KK_ACC_FLOAT>(v[0]);
      ev.v[1] += static_cast<KK_ACC_FLOAT>(v[1]);
      ev.v[2] += static_cast<KK_ACC_FLOAT>(v[2]);
      ev.v[3] += static_cast<KK_ACC_FLOAT>(v[3]);
      ev.v[4] += static_cast<KK_ACC_FLOAT>(v[4]);
      ev.v[5] += static_cast<KK_ACC_FLOAT>(v[5]);
    }

    if (vflag_atom) {
      for (int n = 0; n < 6; n++) {
        const KK_ACC_FLOAT vthird = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[n]);
        a_vatom(i,n) += vthird;
        a_vatom(j,n) += vthird;
        a_vatom(k,n) += vthird;
      }
    }
  }
}

namespace LAMMPS_NS {
template class PairNb3bHarmonicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairNb3bHarmonicKokkos<LMPHostType>;
#endif
}
