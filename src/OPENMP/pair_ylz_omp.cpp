// clang-format off
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

#include "pair_ylz_omp.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairYLZOMP::PairYLZOMP(LAMMPS *lmp) : PairYLZ(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairYLZOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;
  int flag_thr = 0;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag) reduction(+:flag_thr)
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
          eval<1, 1, 1>(ifrom, ito, thr, flag_thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr, flag_thr);
      } else {
        if (force->newton_pair)
          eval<1, 0, 1>(ifrom, ito, thr, flag_thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr, flag_thr);
      }
    } else {
      if (force->newton_pair)
        eval<0, 0, 1>(ifrom, ito, thr, flag_thr);
      else
        eval<0, 0, 0>(ifrom, ito, thr, flag_thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region

  if (flag_thr) {
    int flag_all;
    MPI_Allreduce(&flag_thr, &flag_all, 1, MPI_INT, MPI_MAX, world);
    if (flag_all) error->all(FLERR, "All atoms for pair style ylz must be ellipsoids");
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairYLZOMP::eval(int iifrom, int iito, ThrData *const thr, int &flag)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double evdwl, one_eng, rsq, factor_lj;
  double fforce[3], ttor[3], rtor[3], r12[3];
  double a1[3][3], a2[3][3];
  int *ilist, *jlist, *numneigh, **firstneigh;
  double *iquat, *jquat;

  evdwl = 0.0;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  const int *_noalias const ellipsoid = atom->ellipsoid;
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  auto *_noalias const tor = (dbl3_t *) thr->get_torque()[0];
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *_noalias const special_lj = force->special_lj;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    itype = type[i];

    if (ellipsoid[i] < 0) {
      ++flag;
      continue;
    }

    iquat = bonus[ellipsoid[i]].quat;
    MathExtra::quat_to_mat_trans(iquat, a1);

    jlist = firstneigh[i];
    jnum = numneigh[i];

    double fxtmp, fytmp, fztmp, t1tmp, t2tmp, t3tmp;
    fxtmp = fytmp = fztmp = t1tmp = t2tmp = t3tmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      r12[0] = x[j].x - x[i].x;
      r12[1] = x[j].y - x[i].y;
      r12[2] = x[j].z - x[i].z;
      rsq = MathExtra::dot3(r12, r12);
      jtype = type[j];

      if (ellipsoid[j] < 0) {
        ++flag;
        continue;
      }

      if (rsq < cutsq[itype][jtype]) {
        jquat = bonus[ellipsoid[j]].quat;
        MathExtra::quat_to_mat_trans(jquat, a2);
        one_eng = ylz_analytic(i, j, a1, a2, r12, rsq, fforce, ttor, rtor);

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        fxtmp += fforce[0];
        fytmp += fforce[1];
        fztmp += fforce[2];
        t1tmp += ttor[0];
        t2tmp += ttor[1];
        t3tmp += ttor[2];

        if (NEWTON_PAIR || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j].x -= fforce[0];
          f[j].y -= fforce[1];
          f[j].z -= fforce[2];
          tor[j].x += rtor[0];
          tor[j].y += rtor[1];
          tor[j].z += rtor[2];
        }

        if (EFLAG) evdwl = factor_lj * one_eng;

        if (EVFLAG)
          ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0,
                           fforce[0], fforce[1], fforce[2],
                           -r12[0], -r12[1], -r12[2], thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    tor[i].x += t1tmp;
    tor[i].y += t2tmp;
    tor[i].z += t3tmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairYLZOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairYLZ::memory_usage();
  return bytes;
}
