/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_lepton_coul_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "Lepton.h"
#include "lepton_utils.h"
#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLeptonCoulOMP::PairLeptonCoulOMP(LAMMPS *lmp) : PairLeptonCoul(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLeptonCoulOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair)
          eval<1, 1, 1>(ifrom, ito, thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair)
          eval<1, 0, 1>(ifrom, ito, thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair)
        eval<0, 0, 1>(ifrom, ito, thr);
      else
        eval<0, 0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLeptonCoulOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const auto *_noalias const q = atom->q;
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *_noalias const special_coul = force->special_coul;

  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  double fxtmp, fytmp, fztmp;

  const double q2e = sqrt(force->qqrd2e);

  std::vector<Lepton::CompiledExpression> pairforce;
  std::vector<Lepton::CompiledExpression> pairpot;
  std::vector<std::pair<bool, bool>> have_q;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, Pointers::lmp), functions);
      pairforce.emplace_back(parsed.differentiate("r").createCompiledExpression());
      if (EFLAG) pairpot.emplace_back(parsed.createCompiledExpression());
      pairforce.back().getVariableReference("r");
      have_q.emplace_back(true, true);

      // check if there are references to charges
      try {
        pairforce.back().getVariableReference("qi");
      } catch (std::exception &) {
        have_q.back().first = false;
      }
      try {
        pairforce.back().getVariableReference("qj");
      } catch (std::exception &) {
        have_q.back().second = false;
      }
    }
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      const int jtype = type[j];

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        const int idx = type2expression[itype][jtype];
        pairforce[idx].getVariableReference("r") = r;
        if (have_q[idx].first) pairforce[idx].getVariableReference("qi") = q2e * q[i];
        if (have_q[idx].second) pairforce[idx].getVariableReference("qj") = q2e * q[j];
        const double fpair = -pairforce[idx].evaluate() / r * factor_coul;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        double ecoul = 0.0;
        if (EFLAG) {
          pairpot[idx].getVariableReference("r") = r;
          if (have_q[idx].first) pairpot[idx].getVariableReference("qi") = q2e * q[i];
          if (have_q[idx].second) pairpot[idx].getVariableReference("qj") = q2e * q[j];
          ecoul = pairpot[idx].evaluate();
          ecoul *= factor_coul;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, ecoul, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLeptonCoulOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLeptonCoul::memory_usage();

  return bytes;
}
