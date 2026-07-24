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

#include "pair_lj_cut_coul_long2_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "neigh_list_kokkos.h"

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ----------------------------------------------------------------------
   Hand-written force-only kernel for lj/cut/coul/long, HALFTHREAD + newton.

   Warp-per-atom over a HALF list: one warp's vector lanes split atom i's
   neighbor loop, the i-force is reduced across the lanes, and each pair
   scatters -f to j with an atomic.  This is the configuration the stock
   neigh/thread path cannot express (it requires newton off + a full list,
   doubling the pair count); a half list keeps the pair count halved while
   still getting the memory-level parallelism that hides the x[j] gather
   latency.  Coulomb uses the analytical Ewald correction (no table lookup).
   Parameters travel in the functor (constant memory on CUDA) via the stack
   arrays, as in the STACKPARAMS path of the generic kernel.
------------------------------------------------------------------------- */

namespace {

template<class DeviceType>
struct LJCL2Force {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;

  int nlocal, newton_pair, inum, atoms_per_team;
  KK_FLOAT g_ewald, qqrd2e;
  KK_FLOAT special_lj[4], special_coul[4];
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  // Force contribution of one (i,jj) neighbor pair: folds the i-force into
  // (fx,fy,fz) and atomic-scatters -f to j.  Pulled out of the neighbor loop so
  // the 2-way unroll below can inline it twice; two independent inlined copies
  // let the compiler keep both (irregular) x(j) gathers in flight at once.
  KOKKOS_INLINE_FUNCTION
  void pair_contrib(const int i, const KK_FLOAT xtmp, const KK_FLOAT ytmp,
                    const KK_FLOAT ztmp, const KK_FLOAT qtmp, const int itype,
                    const int jj, KK_ACC_FLOAT& fx, KK_ACC_FLOAT& fy,
                    KK_ACC_FLOAT& fz) const {
    const int jraw = d_neighbors(i,jj);
    const int sb = (jraw >> SBBITS) & 3;
    const KK_FLOAT factor_lj = special_lj[sb];
    const KK_FLOAT factor_coul = special_coul[sb];
    const int j = jraw & NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < m_cutsq[itype][jtype]) {
      // one reciprocal-sqrt supplies rinv, r2inv, and r (was 2 divides + 1 sqrt)
      const KK_FLOAT rinv = Kokkos::rsqrt(rsq);
      const KK_FLOAT r2inv = rinv*rinv;
      KK_FLOAT fpair = static_cast<KK_FLOAT>(0.0);

      if (rsq < m_cut_ljsq[itype][jtype]) {
        const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
        fpair += factor_lj *
          r6inv*(m_params[itype][jtype].lj1*r6inv - m_params[itype][jtype].lj2) * r2inv;
      }

      if (rsq < m_cut_coulsq[itype][jtype]) {
        // analytical Ewald real-space correction (Abramowitz-Stegun 7.1.26)
        const KK_FLOAT r = rsq*rinv;
        const KK_FLOAT grij = g_ewald * r;
        const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
        const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
          (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
        const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
                              t * (static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+
                              t * static_cast<KK_FLOAT>(A5))))) * expm2;
        const KK_FLOAT prefactor = qqrd2e * qtmp*q(j) * rinv;
        KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
        if (factor_coul < static_cast<KK_FLOAT>(1.0))
          forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
        fpair += forcecoul * r2inv;
      }

      fx += static_cast<KK_ACC_FLOAT>(delx*fpair);
      fy += static_cast<KK_ACC_FLOAT>(dely*fpair);
      fz += static_cast<KK_ACC_FLOAT>(delz*fpair);

      if (newton_pair || j < nlocal) {
        Kokkos::atomic_add(&f(j,0), -static_cast<KK_ACC_FLOAT>(delx*fpair));
        Kokkos::atomic_add(&f(j,1), -static_cast<KK_ACC_FLOAT>(dely*fpair));
        Kokkos::atomic_add(&f(j,2), -static_cast<KK_ACC_FLOAT>(delz*fpair));
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
    const int atom_in_team = team.league_rank()*atoms_per_team + team.team_rank();
    if (atom_in_team >= inum) return;

    const int i = d_ilist(atom_in_team);
    const KK_FLOAT xtmp = x(i,0);
    const KK_FLOAT ytmp = x(i,1);
    const KK_FLOAT ztmp = x(i,2);
    const KK_FLOAT qtmp = q(i);
    const int itype = type(i);
    const int jnum = d_numneigh(i);

    KK_ACC_FLOAT fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

    // Vector lanes split the neighbor loop, 2 neighbors per lane per step: the
    // two x(j) gathers of a step are independent, so the compiler overlaps their
    // memory latency (ILP) -- latency hiding that raising occupancy can't reach.
    // Contributions are reduced across lanes into (fxtmp,fytmp,fztmp).
    const int npair = (jnum + 1) >> 1;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, npair),
      [&] (const int k, KK_ACC_FLOAT& fx, KK_ACC_FLOAT& fy, KK_ACC_FLOAT& fz) {
      const int jj0 = k << 1;
      pair_contrib(i, xtmp, ytmp, ztmp, qtmp, itype, jj0, fx, fy, fz);
      const int jj1 = jj0 + 1;
      if (jj1 < jnum)
        pair_contrib(i, xtmp, ytmp, ztmp, qtmp, itype, jj1, fx, fy, fz);
    }, fxtmp, fytmp, fztmp);

    // one lane per atom folds the reduced i-force into the global array
    Kokkos::single(Kokkos::PerThread(team), [&] () {
      Kokkos::atomic_add(&f(i,0), fxtmp);
      Kokkos::atomic_add(&f(i,1), fytmp);
      Kokkos::atomic_add(&f(i,2), fztmp);
    });
  }
};

}  // anonymous namespace

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutCoulLong2Kokkos<DeviceType>::PairLJCutCoulLong2Kokkos(LAMMPS *lmp) :
  PairLJCutCoulLongKokkos<DeviceType>(lmp) {}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  // Optimized path handles only the force-only, HALFTHREAD + newton, stack-param
  // case (the melt hot path).  Everything else -- energy/virial steps for thermo,
  // full or plain-half lists, newton off, many atom types -- uses the proven
  // base-class kernel unchanged.
  if (eflag_in || vflag_in ||
      this->neighflag != HALFTHREAD ||
      !this->force->newton_pair ||
      this->atom->ntypes > MAX_TYPES_STACKPARAMS) {
    PairLJCutCoulLongKokkos<DeviceType>::compute(eflag_in, vflag_in);
    return;
  }

  this->eflag = eflag_in;
  this->vflag = vflag_in;
  this->ev_init(eflag_in, vflag_in, 0);

  this->atomKK->sync(this->execution_space, X_MASK | F_MASK | TYPE_MASK | Q_MASK);
  this->atomKK->modified(this->execution_space, F_MASK);

  LJCL2Force<DeviceType> ff;
  ff.x    = this->atomKK->k_x.template view<DeviceType>();
  ff.f    = this->atomKK->k_f.template view<DeviceType>();
  ff.q    = this->atomKK->k_q.template view<DeviceType>();
  ff.type = this->atomKK->k_type.template view<DeviceType>();

  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
  ff.d_ilist    = k_list->d_ilist;
  ff.d_numneigh = k_list->d_numneigh;
  ff.d_neighbors= k_list->d_neighbors;
  ff.inum       = k_list->inum;

  ff.nlocal      = this->atom->nlocal;
  ff.newton_pair = this->force->newton_pair;
  ff.g_ewald     = static_cast<KK_FLOAT>(this->g_ewald);
  ff.qqrd2e      = static_cast<KK_FLOAT>(this->force->qqrd2e);
  for (int k = 0; k < 4; k++) {
    ff.special_lj[k]   = static_cast<KK_FLOAT>(this->force->special_lj[k]);
    ff.special_coul[k] = static_cast<KK_FLOAT>(this->force->special_coul[k]);
  }

  const int nt = this->atom->ntypes;
  for (int i = 1; i <= nt; i++)
    for (int j = 1; j <= nt; j++) {
      ff.m_params[i][j]      = this->m_params[i][j];
      ff.m_cutsq[i][j]       = this->m_cutsq[i][j];
      ff.m_cut_ljsq[i][j]    = this->m_cut_ljsq[i][j];
      ff.m_cut_coulsq[i][j]  = this->m_cut_coulsq[i][j];
    }

  // warp-per-atom: a full warp's vector lanes split each atom's neighbor loop;
  // atoms_per_team of them share one thread block (4 x 32 = 128-thread blocks).
  // These two constants are the launch-shape tuning knobs -- adjust and rebuild
  // to sweep on new hardware (the stock threads/per/atom and pair/team/size
  // package options are gated behind neigh/thread on and so are unavailable to
  // this half-list, newton-on path).
  const int vector_length = 32;
  const int atoms_per_team = 4;
  ff.atoms_per_team = atoms_per_team;
  const int nteams = (k_list->inum + atoms_per_team - 1) / atoms_per_team;

  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  policy_t policy(nteams, atoms_per_team, vector_length);
  Kokkos::parallel_for("PairLJCutCoulLong2::force", policy, ff);
}

namespace LAMMPS_NS {
template class PairLJCutCoulLong2Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutCoulLong2Kokkos<LMPHostType>;
#endif
}
