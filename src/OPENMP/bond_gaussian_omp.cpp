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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include "bond_gaussian_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 2.3e-308;

/* ---------------------------------------------------------------------- */

BondGaussianOMP::BondGaussianOMP(class LAMMPS *lmp)
  : BondGaussian(lmp), ThrOMP(lmp,THR_BOND)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void BondGaussianOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nbondlist;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
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
          if (force->newton_bond) eval<1,1,1>(ifrom, ito, thr);
          else eval<1,1,0>(ifrom, ito, thr);
        } else {
          if (force->newton_bond) eval<1,0,1>(ifrom, ito, thr);
          else eval<1,0,0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_bond) eval<0,0,1>(ifrom, ito, thr);
        else eval<0,0,0>(ifrom, ito, thr);
      }
    }
    thr->timer(Timer::BOND);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void BondGaussianOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr;
  double prefactor,exponent,g_i,sum_g_i,sum_numerator;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int3_t * _noalias const bondlist = (int3_t *) neighbor->bondlist[0];
  const int nlocal = atom->nlocal;
  ebond = 0.0;

  for (n = nfrom; n < nto; n++) {
    i1 = bondlist[n].a;
    i2 = bondlist[n].b;
    type = bondlist[n].t;

    delx = x[i1].x - x[i2].x;
    dely = x[i1].y - x[i2].y;
    delz = x[i1].z - x[i2].z;

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);

    sum_g_i = 0.0;
    sum_numerator = 0.0;
    for (int i = 0; i < nterms[type]; i++) {
      dr = r - r0[type][i];
      prefactor = (alpha[type][i] / (width[type][i] * sqrt(MY_PI2)));
      exponent = -2.0 * dr * dr / (width[type][i] * width[type][i]);
      g_i = prefactor * exp(exponent);
      sum_g_i += g_i;
      sum_numerator += g_i * dr / (width[type][i] * width[type][i]);
    }

    // avoid overflow
    if (sum_g_i < sum_numerator * SMALL) sum_g_i = sum_numerator * SMALL;

    // force & energy
    if (r > 0.0)
      fbond = -4.0 * (force->boltz * bond_temperature[type]) * (sum_numerator / sum_g_i) / r;
    else
      fbond = 0.0;

    if (EFLAG) ebond = -(force->boltz * bond_temperature[type]) * log(sum_g_i);

    // apply force to each of 2 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += delx*fbond;
      f[i1].y += dely*fbond;
      f[i1].z += delz*fbond;
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x -= delx*fbond;
      f[i2].y -= dely*fbond;
      f[i2].z -= delz*fbond;
    }

    if (EVFLAG) ev_tally_thr(this,i1,i2,nlocal,NEWTON_BOND,
                             ebond,fbond,delx,dely,delz,thr);
  }
}
