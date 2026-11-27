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

#include "angle_cosine_buck6d_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "pair.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

AngleCosineBuck6dOMP::AngleCosineBuck6dOMP(LAMMPS *lmp) :
    AngleCosineBuck6d(lmp), ThrOMP(lmp, THR_ANGLE | THR_CHARMM)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void AngleCosineBuck6dOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // ensure pair->ev_tally() will use 1-3 virial contribution

  if (vflag_global == VIRIAL_FDOTR)
    force->pair->vflag_either = force->pair->vflag_global = 1;

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
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, cvatom, thr);

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
void AngleCosineBuck6dOMP::eval(int nfrom, int nto, ThrData *const thr)
{
  int i1, i2, i3, n, type, itype, jtype;
  double delx1, dely1, delz1, delx2, dely2, delz2;
  double eangle, f1[3], f3[3];
  double rsq1, rsq2, r1, r2, c, s, a, a11, a12, a22;
  double tk;

  // extra lj variables
  double delx3, dely3, delz3, rsq3, r3;
  double rexp, r32inv, r36inv, r314inv, forcebuck6d, fpair;
  double term1, term2, term3, term4, term5, ebuck6d, evdwl;
  double rcu, rqu, sme, smf;

  eangle = 0.0;

  const double *const *const x = atom->x;
  double **const f = thr->get_f();
  const int *const *const anglelist = neighbor->anglelist;
  const int nlocal = atom->nlocal;
  const int newton_pair = force->newton_pair;
  const int *const atomtype = atom->type;

  for (n = nfrom; n < nto; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    r2 = sqrt(rsq2);

    // c = cosine of angle

    c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
    c /= r1 * r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c * c);
    if (s < SMALL) s = SMALL;
    s = 1.0 / s;

    // force & energy

    // explicit lj-contribution

    itype = atomtype[i1];
    jtype = atomtype[i3];

    delx3 = x[i1][0] - x[i3][0];
    dely3 = x[i1][1] - x[i3][1];
    delz3 = x[i1][2] - x[i3][2];
    rsq3 = delx3 * delx3 + dely3 * dely3 + delz3 * delz3;

    r32inv = 0.0;
    if (rsq3 < cut_ljsq[itype][jtype]) {
      r3 = sqrt(rsq3);
      r32inv = 1.0 / rsq3;
      r36inv = r32inv * r32inv * r32inv;
      r314inv = r36inv * r36inv * r32inv;
      rexp = exp(-r3 * buck6d2[itype][jtype]);
      term1 = buck6d3[itype][jtype] * r36inv;
      term2 = buck6d4[itype][jtype] * r314inv;
      term3 = term2 * term2;
      term4 = 1.0 / (1.0 + term2);
      term5 = 1.0 / (1.0 + 2.0 * term2 + term3);
      forcebuck6d = buck6d1[itype][jtype] * buck6d2[itype][jtype] * r3 * rexp;
      forcebuck6d -= term1 * (6.0 * term4 - term5 * 14.0 * term2);
      ebuck6d = buck6d1[itype][jtype] * rexp - term1 * term4;

      // smoothing term
      if (rsq3 > rsmooth_sq[itype][jtype]) {
        rcu = r3 * rsq3;
        rqu = rsq3 * rsq3;
        sme = c5[itype][jtype] * rqu * r3 + c4[itype][jtype] * rqu + c3[itype][jtype] * rcu +
            c2[itype][jtype] * rsq3 + c1[itype][jtype] * r3 + c0[itype][jtype];
        smf = 5.0 * c5[itype][jtype] * rqu + 4.0 * c4[itype][jtype] * rcu +
            3.0 * c3[itype][jtype] * rsq3 + 2.0 * c2[itype][jtype] * r3 + c1[itype][jtype];
        forcebuck6d = forcebuck6d * sme + ebuck6d * smf;
        ebuck6d *= sme;
      }
    } else
      forcebuck6d = 0.0;

    // add forces of additional LJ interaction

    fpair = forcebuck6d * r32inv;
    if (newton_pair || i1 < nlocal) {
      f[i1][0] += delx3 * fpair;
      f[i1][1] += dely3 * fpair;
      f[i1][2] += delz3 * fpair;
    }
    if (newton_pair || i3 < nlocal) {
      f[i3][0] -= delx3 * fpair;
      f[i3][1] -= dely3 * fpair;
      f[i3][2] -= delz3 * fpair;
    }

    evdwl = 0.0;
    if (EFLAG) {
      if (rsq3 < cut_ljsq[itype][jtype]) { evdwl = ebuck6d - offset[itype][jtype]; }
    }

    // update pair energy and virial via per-thread pair accumulators

    if (EVFLAG)
      ev_tally_thr(force->pair, i1, i3, nlocal, newton_pair, evdwl, 0.0, fpair, delx3, dely3, delz3,
                   thr);

    tk = multiplicity[type] * acos(c) - th0[type];

    if (EFLAG) eangle = k[type] * (1.0 + cos(tk));

    a = k[type] * multiplicity[type] * sin(tk) * s;

    a11 = a * c / rsq1;
    a12 = -a / (r1 * r2);
    a22 = a * c / rsq2;

    f1[0] = a11 * delx1 + a12 * delx2;
    f1[1] = a11 * dely1 + a12 * dely2;
    f1[2] = a11 * delz1 + a12 * delz2;
    f3[0] = a22 * delx2 + a12 * delx1;
    f3[1] = a22 * dely2 + a12 * dely1;
    f3[2] = a22 * delz2 + a12 * delz1;

    // apply force to each of 3 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (EVFLAG)
      ev_tally_thr(this, i1, i2, i3, nlocal, NEWTON_BOND, eangle, f1, f3, delx1, dely1, delz1, delx2,
                   dely2, delz2, thr);
  }
}
