// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
------------------------------------------------------------------------- */

#include "pair_body_rounded_polygon_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr int NFNC = 12;    // same as in PairBodyRoundedPolygon

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygonOMP::PairBodyRoundedPolygonOMP(LAMMPS *lmp) :
    PairBodyRoundedPolygon(lmp), ThrOMP(lmp, THR_PAIR), fnc_thr(nullptr), nmax_thr(0)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygonOMP::~PairBodyRoundedPolygonOMP()
{
  memory->destroy(fnc_thr);
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygonOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow the per-atom lists if necessary and initialize

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum, nmax, "pair:dnum");
    memory->create(dfirst, nmax, "pair:dfirst");
    memory->create(ednum, nmax, "pair:ednum");
    memory->create(edfirst, nmax, "pair:edfirst");
    memory->create(enclosing_radius, nmax, "pair:enclosing_radius");
    memory->create(rounded_radius, nmax, "pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (int i = 0; i < nall; i++) dnum[i] = ednum[i] = 0;

  // body2space modifies the shared vertex and edge lists,
  // so it is called for all bodies before the parallel region

  int *body = atom->body;
  for (int i = 0; i < nall; i++)
    if (body[i] >= 0) body2space(i);

  // per-thread forces and torques that do not derive from the energy,
  // and per-thread scratch space

  if (atom->nmax > nmax_fnc) {
    memory->destroy(fnc);
    nmax_fnc = atom->nmax;
    memory->create(fnc, nmax_fnc, NFNC, "pair:fnc");
  }
  if ((atom->nmax > nmax_thr) || ((int) scratch_thr.size() != nthreads)) {
    memory->destroy(fnc_thr);
    nmax_thr = atom->nmax;
    memory->create(fnc_thr, nthreads * nmax_thr, NFNC, "pair:fnc_thr");
    scratch_thr.resize(nthreads);
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    double **fnc_t = fnc_thr + tid * nmax_thr;
    for (int i = 0; i < nall; i++)
      for (int k = 0; k < NFNC; k++) fnc_t[i][k] = 0.0;

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1, 1, 1>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
        else eval<1, 1, 0>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
      } else {
        if (force->newton_pair) eval<1, 0, 1>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
        else eval<1, 0, 0>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
      }
    } else {
      eval<0, 0, 0>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region

  // sum the per-thread forces and torques that do not derive from the energy

  for (int i = 0; i < nall; i++) {
    for (int k = 0; k < NFNC; k++) {
      double sum = 0.0;
      for (int t = 0; t < nthreads; t++) sum += fnc_thr[t * nmax_thr + i][k];
      fnc[i][k] = sum;
    }
  }

  work_nonconservative();
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBodyRoundedPolygonOMP::eval(int iifrom, int iito, ThrData *const thr,
                                     double **fnc_t, Scratch &s)
{
  double **x = atom->x;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double **f = thr->get_f();
  double **torque = thr->get_torque();
  double *radius = atom->radius;
  int *body = atom->body;
  const int nlocal = atom->nlocal;

  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  int **const firstneigh = list->firstneigh;

  double facc[3];

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    if (body[i] < 0) continue;
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const double radi = radius[i];
    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj] & NEIGHMASK;
      if (body[j] < 0) continue;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;

      // no interaction

      if (sqrt(rsq) > radi + radius[j] + cut_inner) continue;

      double evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0.0;
      pair_interaction(i, j, delx, dely, delz, rsq, x, v, angmom, f, torque, fnc_t, s, evdwl,
                       facc);

      if (EVFLAG)
        ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR, EFLAG ? evdwl : 0.0, 0.0, facc[0],
                         facc[1], facc[2], delx, dely, delz, thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBodyRoundedPolygonOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBodyRoundedPolygon::memory_usage();
  bytes += (double) nmax_thr * scratch_thr.size() * NFNC * sizeof(double);
  return bytes;
}
