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

#include "angle_lepton_omp.h"
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

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

AngleLeptonOMP::AngleLeptonOMP(class LAMMPS *lmp) : AngleLepton(lmp), ThrOMP(lmp, THR_ANGLE)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void AngleLeptonOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nanglelist;

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
void AngleLeptonOMP::eval(int nfrom, int nto, ThrData *const thr)
{
  std::vector<Lepton::CompiledExpression> angleforce;
  std::vector<Lepton::CompiledExpression> anglepot;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, Pointers::lmp));
      angleforce.emplace_back(parsed.differentiate("theta").createCompiledExpression());
      if (EFLAG) anglepot.emplace_back(parsed.createCompiledExpression());
    }
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int4_t *_noalias const anglelist = (int4_t *) neighbor->anglelist[0];
  const int nlocal = atom->nlocal;

  for (int n = nfrom; n < nto; n++) {
    const int i1 = anglelist[n].a;
    const int i2 = anglelist[n].b;
    const int i3 = anglelist[n].c;
    const int type = anglelist[n].t;

    // 1st bond

    const double delx1 = x[i1].x - x[i2].x;
    const double dely1 = x[i1].y - x[i2].y;
    const double delz1 = x[i1].z - x[i2].z;

    const double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    const double r1 = sqrt(rsq1);

    // 2nd bond

    const double delx2 = x[i3].x - x[i2].x;
    const double dely2 = x[i3].y - x[i2].y;
    const double delz2 = x[i3].z - x[i2].z;

    const double rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    const double r2 = sqrt(rsq2);

    // angle (cos and sin)

    double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
    c /= r1 * r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    double s = sqrt(1.0 - c * c);
    if (s < SMALL) s = SMALL;
    s = 1.0 / s;

    // force and energy

    const double dtheta = acos(c) - theta0[type];
    const int idx = type2expression[type];
    angleforce[idx].getVariableReference("theta") = dtheta;

    const double a = -angleforce[idx].evaluate() * s;
    const double a11 = a * c / rsq1;
    const double a12 = -a / (r1 * r2);
    const double a22 = a * c / rsq2;

    double f1[3], f3[3];
    f1[0] = a11 * delx1 + a12 * delx2;
    f1[1] = a11 * dely1 + a12 * dely2;
    f1[2] = a11 * delz1 + a12 * delz2;
    f3[0] = a22 * delx2 + a12 * delx1;
    f3[1] = a22 * dely2 + a12 * dely1;
    f3[2] = a22 * delz2 + a12 * delz1;

    // apply force to each of 3 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += f1[0];
      f[i1].y += f1[1];
      f[i1].z += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x -= f1[0] + f3[0];
      f[i2].y -= f1[1] + f3[1];
      f[i2].z -= f1[2] + f3[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3].x += f3[0];
      f[i3].y += f3[1];
      f[i3].z += f3[2];
    }

    double eangle = 0.0;
    if (EFLAG) {
      anglepot[idx].getVariableReference("theta") = dtheta;
      eangle = anglepot[idx].evaluate() - offset[type];
    }
    if (EVFLAG)
      ev_tally_thr(this, i1, i2, i3, nlocal, NEWTON_BOND, eangle, f1, f3, delx1, dely1, delz1,
                   delx2, dely2, delz2, thr);
  }
}
