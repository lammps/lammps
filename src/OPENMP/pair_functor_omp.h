/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_FUNCTOR_OMP_H
#define LMP_PAIR_FUNCTOR_OMP_H

// OpenMP-threaded driver for the FUNCTOR pair styles.  It is the ThrOMP analog
// of PairFunctor<EVAL,COUL>: it reuses the same evaluator EVAL and Coulomb
// policy COUL and the same per-pair parameter storage from the base class, but
// the neighbor loop is split across threads, each thread accumulates into its
// own force slice (ThrData::get_f), and energy/virial are tallied through
// ev_tally_thr and combined by reduce_thr.  A concrete style is a thin subclass,
// e.g. PairLJCutFunctorOMP : public PairFunctorOMP<EvaluatorLJCut, CoulNone>.

#include "pair_functor.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "thr_omp.h"

#include "functor_vec.h"

#include "omp_compat.h"

#include <cstddef>

namespace LAMMPS_NS {

template <class EVAL, class COUL>
class PairFunctorOMP : public PairFunctor<EVAL, COUL>, public ThrOMP {
 protected:
  using Base = PairFunctor<EVAL, COUL>;
  using Param = typename EVAL::Param;

 public:
  PairFunctorOMP(LAMMPS *lmp) : Base(lmp), ThrOMP(lmp, THR_PAIR)
  {
    this->suffix_flag |= Suffix::OMP;
    this->respa_enable = 0;
  }

  void compute(int eflag, int vflag) override
  {
    this->ev_init(eflag, vflag);

    const int nall = this->atom->nlocal + this->atom->nghost;
    const int nthreads = this->comm->nthreads;
    const int inum = this->list->inum;
    const int newton = this->force->newton_pair;
    const int ctable = (COUL::has_table && this->ncoultablebits != 0) ? 1 : 0;
    const int special = (this->neighbor->special_flag[1] == 2 ||
                         this->neighbor->special_flag[2] == 2 ||
                         this->neighbor->special_flag[3] == 2)
        ? 1
        : 0;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
    {
      int ifrom, ito, tid;
      loop_setup_thr(ifrom, ito, tid, inum, nthreads);
      ThrData *thr = fix->get_thr(tid);
      thr->timer(Timer::START);
      ev_setup_thr(eflag, vflag, nall, this->eatom, this->vatom, nullptr, thr);

      if (this->evflag) {
        if (eflag) {
          if (newton) dispatch<1, 1, 1>(ctable, special, ifrom, ito, thr);
          else dispatch<1, 1, 0>(ctable, special, ifrom, ito, thr);
        } else {
          if (newton) dispatch<1, 0, 1>(ctable, special, ifrom, ito, thr);
          else dispatch<1, 0, 0>(ctable, special, ifrom, ito, thr);
        }
      } else {
        if (newton) dispatch<0, 0, 1>(ctable, special, ifrom, ito, thr);
        else dispatch<0, 0, 0>(ctable, special, ifrom, ito, thr);
      }

      thr->timer(Timer::PAIR);
      reduce_thr(this, eflag, vflag, thr);
    }
  }

  double memory_usage() override { return memory_usage_thr() + Base::memory_usage(); }

 protected:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void dispatch(int ctable, int special, int ifrom, int ito, ThrData *thr)
  {
    if (ctable) {
      if (special) eval<EVFLAG, EFLAG, NEWTON_PAIR, 1, 1>(ifrom, ito, thr);
      else eval<EVFLAG, EFLAG, NEWTON_PAIR, 1, 0>(ifrom, ito, thr);
    } else {
      if (special) eval<EVFLAG, EFLAG, NEWTON_PAIR, 0, 1>(ifrom, ito, thr);
      else eval<EVFLAG, EFLAG, NEWTON_PAIR, 0, 0>(ifrom, ito, thr);
    }
  }

  template <int EVFLAG, int EFLAG, int NEWTON_PAIR, int CTABLE, int SPECIAL>
  void eval(int iifrom, int iito, ThrData *const thr)
  {
    // local aliases for members inherited through the dependent template base
    const auto *_noalias xx = (functor::vec3_t *) this->atom->x[0];
    auto *_noalias ff = (functor::vec3_t *) thr->get_f()[0];
    const int *_noalias type = this->atom->type;
    const int nlocal = this->atom->nlocal;
    const double *_noalias special_lj = this->force->special_lj;
    const Param *_noalias params = this->params;
    const int np = this->nparams;
    const auto &coul = this->coul;

    [[maybe_unused]] const double *_noalias q = nullptr;
    [[maybe_unused]] const double *_noalias special_coul = nullptr;
    [[maybe_unused]] const double *_noalias coul_cutsq = nullptr;
    if constexpr (COUL::needs_charge) {
      q = this->atom->q;
      special_coul = this->force->special_coul;
    }
    if constexpr (COUL::has_coul) coul_cutsq = coul.cut_coulsq;

    const int *_noalias ilist = this->list->ilist;
    const int *_noalias numneigh = this->list->numneigh;
    const int *const *_noalias firstneigh = this->list->firstneigh;

    for (int ii = iifrom; ii < iito; ++ii) {
      const int i = ilist[ii];
      const double xtmp = xx[i].x;
      const double ytmp = xx[i].y;
      const double ztmp = xx[i].z;
      const int itype = type[i];
      const Param *_noalias pi = params + (std::size_t) itype * np;
      [[maybe_unused]] const double *_noalias coul_cutsq_row = nullptr;
      if constexpr (COUL::has_coul) coul_cutsq_row = coul_cutsq + (std::size_t) itype * np;
      [[maybe_unused]] double qi = 0.0;
      if constexpr (COUL::needs_charge) qi = q[i];

      const int *_noalias jlist = firstneigh[i];
      const int jnum = numneigh[i];

      double tmpfx = 0.0, tmpfy = 0.0, tmpfz = 0.0;

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        double factor_lj = 1.0;
        [[maybe_unused]] double factor_coul = 1.0;
        if constexpr (SPECIAL) {
          const int sb = this->sbmask(j);
          factor_lj = special_lj[sb];
          if constexpr (COUL::has_coul) factor_coul = special_coul[sb];
        }
        j &= NEIGHMASK;

        const double delx = xtmp - xx[j].x;
        const double dely = ytmp - xx[j].y;
        const double delz = ztmp - xx[j].z;
        double rsq = delx * delx + dely * dely + delz * delz;
        // core-shell policies bump rsq by a tiny epsilon (compile-time 0.0 for
        // every non-cs policy, so this folds away); see PairFunctor::eval
        if constexpr (COUL::rsq_epsilon != 0.0) rsq += COUL::rsq_epsilon;
        const int jtype = type[j];
        const Param &p = pi[jtype];

        double fpair = 0.0, evdwl = 0.0, ecoul = 0.0;
        bool within;

        if constexpr (!COUL::has_coul) {
          within = (rsq < p.cutsq);
          if (within) {
            const auto v = EVAL::template pair<EFLAG>(rsq, p, factor_lj);
            fpair = v.fpair;
            evdwl = v.energy;
          }
        } else {
          const double cut_coulsq = coul_cutsq_row[jtype];
          const double outer = (p.cutsq > cut_coulsq) ? p.cutsq : cut_coulsq;
          within = (rsq < outer);
          if (within) {
            if (rsq < p.cutsq) {
              const auto v = EVAL::template pair<EFLAG>(rsq, p, factor_lj);
              fpair = v.fpair;
              evdwl = v.energy;
            }
            if (rsq < cut_coulsq) {
              const auto c = coul.template eval_coul<EFLAG, CTABLE>(rsq, qi, q[j], factor_coul);
              fpair += c.fpair;
              ecoul = c.energy;
            }
          }
        }

        if (within) {
          tmpfx += delx * fpair;
          tmpfy += dely * fpair;
          tmpfz += delz * fpair;
          if (NEWTON_PAIR || j < nlocal) {
            ff[j].x -= delx * fpair;
            ff[j].y -= dely * fpair;
            ff[j].z -= delz * fpair;
          }
          if (EVFLAG)
            ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, ecoul, fpair, delx, dely, delz,
                         thr);
        }
      }

      ff[i].x += tmpfx;
      ff[i].y += tmpfy;
      ff[i].z += tmpfz;
    }
  }
};

}    // namespace LAMMPS_NS

#endif
