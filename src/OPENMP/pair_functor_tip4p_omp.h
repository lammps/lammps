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

#ifndef LMP_PAIR_FUNCTOR_TIP4P_OMP_H
#define LMP_PAIR_FUNCTOR_TIP4P_OMP_H

// OpenMP-threaded driver for the TIP4P FUNCTOR pair styles.  It is the ThrOMP
// analog of PairFunctorTIP4P<EVAL>: it reuses the same evaluator EVAL, the
// CoulLong member, the per-pair parameter store, and the M-site geometry of the
// serial base, but splits the neighbor loop across threads, accumulates into a
// per-thread force slice (ThrData::get_f), and tallies energy/virial through
// ev_tally_thr / ev_tally_list_thr.  The M-site caches (hneigh/newsite, owned by
// the base) are filled with a thread-safe update order: the "computed" sentinel
// hneigh[i][0] is written last, so a racing thread re-derives the site rather
// than reading a half-written entry (worst case is a repeated calculation with
// identical result, never a race).  A concrete style is a thin subclass, e.g.
// PairLJCutTIP4PLongFunctorOMP : public PairFunctorTIP4POMP<EvaluatorLJCut>.

#include "pair_functor_tip4p.h"

#include "ewald_const.h"
#include "suffix.h"
#include "thr_omp.h"

#include "omp_compat.h"

#include <cmath>
#include <cstddef>

namespace LAMMPS_NS {

template <class EVAL> class PairFunctorTIP4POMP : public PairFunctorTIP4P<EVAL>, public ThrOMP {
 protected:
  using Base = PairFunctorTIP4P<EVAL>;
  using Param = typename EVAL::Param;

 public:
  PairFunctorTIP4POMP(LAMMPS *lmp) : Base(lmp), ThrOMP(lmp, THR_PAIR)
  {
    this->suffix_flag |= Suffix::OMP;
    this->respa_enable = 0;
    // no_virial_fdotr_compute is set by the base ctor (M-site virial is explicit)
  }

  void compute(int eflag, int vflag) override
  {
    this->ev_init(eflag, vflag);

    const int nlocal = this->atom->nlocal;
    const int nall = nlocal + this->atom->nghost;

    // (re)allocate the M-site caches owned by the base, and reset the validity
    // markers (hneigh[i][0] = -1 after reneighboring, hneigh[i][2] = 0 each step)
    if (this->atom->nmax > this->nmax) {
      this->nmax = this->atom->nmax;
      this->memory->destroy(this->hneigh);
      this->memory->create(this->hneigh, this->nmax, 3, "pair:hneigh");
      this->memory->destroy(this->newsite);
      this->memory->create(this->newsite, this->nmax, 3, "pair:newsite");
    }
    if (this->neighbor->ago == 0)
      for (int i = 0; i < nall; i++) this->hneigh[i][0] = -1;
    for (int i = 0; i < nall; i++) this->hneigh[i][2] = 0;

    const int nthreads = this->comm->nthreads;
    const int inum = this->list->inum;
    const int ctable = (this->coul.ncoultablebits != 0) ? 1 : 0;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag, ctable)
#endif
    {
      int ifrom, ito, tid;
      loop_setup_thr(ifrom, ito, tid, inum, nthreads);
      ThrData *thr = fix->get_thr(tid);
      thr->timer(Timer::START);
      ev_setup_thr(eflag, vflag, nall, this->eatom, this->vatom, nullptr, thr);

      if (ctable) {
        if (this->evflag) {
          if (eflag) {
            if (vflag) eval<1, 1, 1, 1>(ifrom, ito, thr);
            else eval<1, 1, 1, 0>(ifrom, ito, thr);
          } else {
            if (vflag) eval<1, 1, 0, 1>(ifrom, ito, thr);
            else eval<1, 1, 0, 0>(ifrom, ito, thr);
          }
        } else {
          eval<1, 0, 0, 0>(ifrom, ito, thr);
        }
      } else {
        if (this->evflag) {
          if (eflag) {
            if (vflag) eval<0, 1, 1, 1>(ifrom, ito, thr);
            else eval<0, 1, 1, 0>(ifrom, ito, thr);
          } else {
            if (vflag) eval<0, 1, 0, 1>(ifrom, ito, thr);
            else eval<0, 1, 0, 0>(ifrom, ito, thr);
          }
        } else {
          eval<0, 0, 0, 0>(ifrom, ito, thr);
        }
      }

      thr->timer(Timer::PAIR);
      reduce_thr(this, eflag, vflag, thr);
    }
  }

  double memory_usage() override { return memory_usage_thr() + Base::memory_usage(); }

 protected:
  template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG>
  void eval(int iifrom, int iito, ThrData *const thr)
  {
    using namespace EwaldConst;

    double **x = this->atom->x;
    double **f = thr->get_f();
    const double *q = this->atom->q;
    const tagint *tag = this->atom->tag;
    const int *type = this->atom->type;
    const int nlocal = this->atom->nlocal;
    const double *special_coul = this->force->special_coul;
    const double *special_lj = this->force->special_lj;

    const Param *const params = this->params;
    const int np = this->nparams;
    const int typeO = this->typeO;
    const int typeH = this->typeH;
    const double alpha = this->alpha;
    const double qdist = this->qdist;
    int **const hneigh = this->hneigh;
    double **const newsite = this->newsite;

    const auto &coul = this->coul;
    const double qqrd2e = coul.qqrd2e;
    const double g_ewald = coul.g_ewald;
    const double tabinnersq = coul.tabinnersq;
    const double cut_coul = coul.cut_coul_global;
    const double cut_coulsq = cut_coul * cut_coul;
    const double cut_coulsqplus = (cut_coul + 2.0 * qdist) * (cut_coul + 2.0 * qdist);

    const int *ilist = this->list->ilist;
    const int *numneigh = this->list->numneigh;
    const int *const *firstneigh = this->list->firstneigh;

    double v[6], fd[3], fO[3], fH[3];
    int vlist[6];

    for (int ii = iifrom; ii < iito; ++ii) {
      const int i = ilist[ii];
      const double qtmp = q[i];
      const double xtmp = x[i][0], ytmp = x[i][1], ztmp = x[i][2];
      const int itype = type[i];
      const Param *pi = params + (std::size_t) itype * np;

      // if atom i is a water O, x1 is the M site; else the atom position.
      // sentinel hneigh[i][0] is written LAST so a racing thread re-derives.
      int iH1 = -1, iH2 = -1;
      double *x1;
      if (itype == typeO) {
        if (hneigh[i][0] < 0) {
          iH1 = this->atom->map(tag[i] + 1);
          iH2 = this->atom->map(tag[i] + 2);
          if (iH1 == -1 || iH2 == -1) this->error->one(FLERR, "TIP4P hydrogen is missing");
          if (type[iH1] != typeH || type[iH2] != typeH)
            this->error->one(FLERR, "TIP4P hydrogen has incorrect atom type");
          iH1 = this->domain->closest_image(i, iH1);
          iH2 = this->domain->closest_image(i, iH2);
          this->compute_newsite(x[i], x[iH1], x[iH2], newsite[i]);
          hneigh[i][2] = 1;
          hneigh[i][1] = iH2;
          hneigh[i][0] = iH1;
        } else {
          iH1 = hneigh[i][0];
          iH2 = hneigh[i][1];
          if (hneigh[i][2] == 0) {
            this->compute_newsite(x[i], x[iH1], x[iH2], newsite[i]);
            hneigh[i][2] = 1;
          }
        }
        x1 = newsite[i];
      } else
        x1 = x[i];

      const int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        const double factor_lj = special_lj[this->sbmask(j)];
        const double factor_coul = special_coul[this->sbmask(j)];
        j &= NEIGHMASK;

        double delx = xtmp - x[j][0];
        double dely = ytmp - x[j][1];
        double delz = ztmp - x[j][2];
        double rsq = delx * delx + dely * dely + delz * delz;
        const int jtype = type[j];
        const Param &p = pi[jtype];

        // LJ interaction based on the true atom separation

        if (rsq < p.cutsq) {
          const auto lj = EVAL::template pair<EFLAG>(rsq, p, factor_lj);
          const double forcelj = lj.fpair;

          fxtmp += delx * forcelj;
          fytmp += dely * forcelj;
          fztmp += delz * forcelj;
          f[j][0] -= delx * forcelj;
          f[j][1] -= dely * forcelj;
          f[j][2] -= delz * forcelj;

          const double evdwl = EFLAG ? lj.energy : 0.0;
          if (EVFLAG)
            ev_tally_thr(this, i, j, nlocal, /* newton_pair = */ 1, evdwl, 0.0, forcelj, delx, dely,
                         delz, thr);
        }

        // Coulomb interaction based on the (possibly M-site-shifted) separation

        if (rsq < cut_coulsqplus) {
          int jH1 = -1, jH2 = -1;
          if (itype == typeO || jtype == typeO) {
            double *x2;
            if (jtype == typeO) {
              if (hneigh[j][0] < 0) {
                jH1 = this->atom->map(tag[j] + 1);
                jH2 = this->atom->map(tag[j] + 2);
                if (jH1 == -1 || jH2 == -1) this->error->one(FLERR, "TIP4P hydrogen is missing");
                if (type[jH1] != typeH || type[jH2] != typeH)
                  this->error->one(FLERR, "TIP4P hydrogen has incorrect atom type");
                jH1 = this->domain->closest_image(j, jH1);
                jH2 = this->domain->closest_image(j, jH2);
                this->compute_newsite(x[j], x[jH1], x[jH2], newsite[j]);
                hneigh[j][2] = 1;
                hneigh[j][1] = jH2;
                hneigh[j][0] = jH1;
              } else {
                jH1 = hneigh[j][0];
                jH2 = hneigh[j][1];
                if (hneigh[j][2] == 0) {
                  this->compute_newsite(x[j], x[jH1], x[jH2], newsite[j]);
                  hneigh[j][2] = 1;
                }
              }
              x2 = newsite[j];
            } else
              x2 = x[j];

            delx = x1[0] - x2[0];
            dely = x1[1] - x2[1];
            delz = x1[2] - x2[2];
            rsq = delx * delx + dely * dely + delz * delz;
          }

          if (rsq < cut_coulsq) {
            const double r2inv = 1.0 / rsq;
            double forcecoul, prefactor = 0.0, erfc = 0.0;
            double fraction = 0.0;
            int itable = 0;
            if (!CTABLE || rsq <= tabinnersq) {
              const double r = sqrt(rsq);
              const double grij = g_ewald * r;
              const double expm2 = exp(-grij * grij);
              const double t = 1.0 / (1.0 + EWALD_P * grij);
              erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
              prefactor = qqrd2e * qtmp * q[j] / r;
              forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
              if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
            } else {
              typename Pair::union_int_float_t rsq_lookup;
              rsq_lookup.f = rsq;
              itable = rsq_lookup.i & coul.ncoulmask;
              itable >>= coul.ncoulshiftbits;
              fraction = ((double) rsq_lookup.f - coul.rtable[itable]) * coul.drtable[itable];
              double table = coul.ftable[itable] + fraction * coul.dftable[itable];
              forcecoul = qtmp * q[j] * table;
              if (factor_coul < 1.0) {
                table = coul.ctable[itable] + fraction * coul.dctable[itable];
                prefactor = qtmp * q[j] * table;
                forcecoul -= (1.0 - factor_coul) * prefactor;
              }
            }

            const double cforce = forcecoul * r2inv;

            // apply Coulomb force; an O atom's force acts on the M site and is
            // partitioned over O and its two H atoms (Feenstra, J Comp Chem 20,
            // 786 (1999)), preserving the total force and torque

            int n = 0, key = 0;

            if (itype != typeO) {
              fxtmp += delx * cforce;
              fytmp += dely * cforce;
              fztmp += delz * cforce;
              if (VFLAG) {
                v[0] = x[i][0] * delx * cforce;
                v[1] = x[i][1] * dely * cforce;
                v[2] = x[i][2] * delz * cforce;
                v[3] = x[i][0] * dely * cforce;
                v[4] = x[i][0] * delz * cforce;
                v[5] = x[i][1] * delz * cforce;
                vlist[n++] = i;
              }
            } else {
              key++;
              fd[0] = delx * cforce;
              fd[1] = dely * cforce;
              fd[2] = delz * cforce;
              fO[0] = fd[0] * (1 - alpha);
              fO[1] = fd[1] * (1 - alpha);
              fO[2] = fd[2] * (1 - alpha);
              fH[0] = 0.5 * alpha * fd[0];
              fH[1] = 0.5 * alpha * fd[1];
              fH[2] = 0.5 * alpha * fd[2];
              fxtmp += fO[0];
              fytmp += fO[1];
              fztmp += fO[2];
              f[iH1][0] += fH[0];
              f[iH1][1] += fH[1];
              f[iH1][2] += fH[2];
              f[iH2][0] += fH[0];
              f[iH2][1] += fH[1];
              f[iH2][2] += fH[2];
              if (VFLAG) {
                const double *xH1 = x[iH1];
                const double *xH2 = x[iH2];
                v[0] = x[i][0] * fO[0] + xH1[0] * fH[0] + xH2[0] * fH[0];
                v[1] = x[i][1] * fO[1] + xH1[1] * fH[1] + xH2[1] * fH[1];
                v[2] = x[i][2] * fO[2] + xH1[2] * fH[2] + xH2[2] * fH[2];
                v[3] = x[i][0] * fO[1] + xH1[0] * fH[1] + xH2[0] * fH[1];
                v[4] = x[i][0] * fO[2] + xH1[0] * fH[2] + xH2[0] * fH[2];
                v[5] = x[i][1] * fO[2] + xH1[1] * fH[2] + xH2[1] * fH[2];
                vlist[n++] = i;
                vlist[n++] = iH1;
                vlist[n++] = iH2;
              }
            }

            if (jtype != typeO) {
              f[j][0] -= delx * cforce;
              f[j][1] -= dely * cforce;
              f[j][2] -= delz * cforce;
              if (VFLAG) {
                v[0] -= x[j][0] * delx * cforce;
                v[1] -= x[j][1] * dely * cforce;
                v[2] -= x[j][2] * delz * cforce;
                v[3] -= x[j][0] * dely * cforce;
                v[4] -= x[j][0] * delz * cforce;
                v[5] -= x[j][1] * delz * cforce;
                vlist[n++] = j;
              }
            } else {
              key += 2;
              fd[0] = -delx * cforce;
              fd[1] = -dely * cforce;
              fd[2] = -delz * cforce;
              fO[0] = fd[0] * (1 - alpha);
              fO[1] = fd[1] * (1 - alpha);
              fO[2] = fd[2] * (1 - alpha);
              fH[0] = 0.5 * alpha * fd[0];
              fH[1] = 0.5 * alpha * fd[1];
              fH[2] = 0.5 * alpha * fd[2];
              f[j][0] += fO[0];
              f[j][1] += fO[1];
              f[j][2] += fO[2];
              f[jH1][0] += fH[0];
              f[jH1][1] += fH[1];
              f[jH1][2] += fH[2];
              f[jH2][0] += fH[0];
              f[jH2][1] += fH[1];
              f[jH2][2] += fH[2];
              if (VFLAG) {
                const double *xH1 = x[jH1];
                const double *xH2 = x[jH2];
                v[0] += x[j][0] * fO[0] + xH1[0] * fH[0] + xH2[0] * fH[0];
                v[1] += x[j][1] * fO[1] + xH1[1] * fH[1] + xH2[1] * fH[1];
                v[2] += x[j][2] * fO[2] + xH1[2] * fH[2] + xH2[2] * fH[2];
                v[3] += x[j][0] * fO[1] + xH1[0] * fH[1] + xH2[0] * fH[1];
                v[4] += x[j][0] * fO[2] + xH1[0] * fH[2] + xH2[0] * fH[2];
                v[5] += x[j][1] * fO[2] + xH1[1] * fH[2] + xH2[1] * fH[2];
                vlist[n++] = j;
                vlist[n++] = jH1;
                vlist[n++] = jH2;
              }
            }

            double ecoul = 0.0;
            if (EFLAG) {
              if (!CTABLE || rsq <= tabinnersq) {
                ecoul = prefactor * erfc;
              } else {
                const double table = coul.etable[itable] + fraction * coul.detable[itable];
                ecoul = qtmp * q[j] * table;
              }
              if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
            }

            if (EVFLAG) ev_tally_list_thr(this, key, vlist, v, ecoul, alpha, thr);
          }
        }
      }

      f[i][0] += fxtmp;
      f[i][1] += fytmp;
      f[i][2] += fztmp;
    }
  }
};

}    // namespace LAMMPS_NS

#endif
