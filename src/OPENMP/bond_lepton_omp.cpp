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
------------------------------------------------------------------------- */

#include "bond_lepton_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

#include "Lepton.h"
#include "lepton_utils.h"
#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondLeptonOMP::BondLeptonOMP(class LAMMPS *_lmp) : BondLepton(_lmp), ThrOMP(_lmp, THR_BOND)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void BondLeptonOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nbondlist;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (inum > 0) {
      if (evflag) {
        if (eflag) {
          if (force->newton_bond)
            eval<1, 1, 1>(ifrom, ito, thr);
          else
            eval<1, 1, 0>(ifrom, ito, thr);
        } else {
          if (force->newton_bond)
            eval<1, 0, 1>(ifrom, ito, thr);
          else
            eval<1, 0, 0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_bond)
          eval<0, 0, 1>(ifrom, ito, thr);
        else
          eval<0, 0, 0>(ifrom, ito, thr);
      }
    }
    thr->timer(Timer::BOND);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void BondLeptonOMP::eval(int nfrom, int nto, ThrData *const thr)
{
  std::vector<Lepton::CompiledExpression> bondforce;
  std::vector<Lepton::CompiledExpression> bondpot;
  std::vector<bool> has_ref;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, Pointers::lmp));
      bondforce.emplace_back(parsed.differentiate("r").createCompiledExpression());
      has_ref.push_back(true);
      try {
        bondforce.back().getVariableReference("r");
      } catch (Lepton::Exception &) {
        has_ref.back() = false;
      }
      if (EFLAG) bondpot.emplace_back(parsed.createCompiledExpression());
    }
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int3_t *_noalias const bondlist = (int3_t *) neighbor->bondlist[0];
  const int nlocal = atom->nlocal;

  for (int n = nfrom; n < nto; n++) {
    const int i1 = bondlist[n].a;
    const int i2 = bondlist[n].b;
    const int type = bondlist[n].t;

    const double delx = x[i1].x - x[i2].x;
    const double dely = x[i1].y - x[i2].y;
    const double delz = x[i1].z - x[i2].z;

    const double rsq = delx * delx + dely * dely + delz * delz;
    const double r = sqrt(rsq);
    const double dr = r - r0[type];
    const int idx = type2expression[type];

    // force and energy

    double fbond = 0.0;
    if (r > 0.0) {
      if (has_ref[idx]) bondforce[idx].getVariableReference("r") = dr;
      fbond = -bondforce[idx].evaluate() / r;
    }

    // apply force to each of 2 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += delx * fbond;
      f[i1].y += dely * fbond;
      f[i1].z += delz * fbond;
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x -= delx * fbond;
      f[i2].y -= dely * fbond;
      f[i2].z -= delz * fbond;
    }

    double ebond = 0.0;
    if (EFLAG) {
      try {
        bondpot[idx].getVariableReference("r") = dr;
      } catch (Lepton::Exception &) {
        ;    // ignore -> constant potential
      }
      ebond = bondpot[idx].evaluate() - offset[type];
    }
    if (EVFLAG)
      ev_tally_thr(this, i1, i2, nlocal, NEWTON_BOND, ebond, fbond, delx, dely, delz, thr);
  }
}
