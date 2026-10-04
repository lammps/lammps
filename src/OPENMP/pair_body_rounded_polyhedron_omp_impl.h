/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)

   Member functions of PairBodyRoundedPolyhedronOMPT, the threaded compute
   of pair style body/rounded/polyhedron and of the pair styles derived from
   it.  Included only by the source files of the threaded styles, which
   instantiate the template for their base style.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_BODY_ROUNDED_POLYHEDRON_OMP_IMPL_H
#define LMP_PAIR_BODY_ROUNDED_POLYHEDRON_OMP_IMPL_H

#include "pair_body_rounded_polyhedron_omp.h"

#include "atom.h"
#include "comm.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"

#include <cmath>

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

template <class BASE>
PairBodyRoundedPolyhedronOMPT<BASE>::PairBodyRoundedPolyhedronOMPT(LAMMPS *lmp) :
    BASE(lmp), ThrOMP(lmp, THR_PAIR), fnc_thr(nullptr), nmax_thr(0)
{
  this->suffix_flag |= Suffix::OMP;
  this->respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

template <class BASE> PairBodyRoundedPolyhedronOMPT<BASE>::~PairBodyRoundedPolyhedronOMPT()
{
  this->memory->destroy(fnc_thr);
}

/* ---------------------------------------------------------------------- */

template <class BASE> void PairBodyRoundedPolyhedronOMPT<BASE>::compute(int eflag, int vflag)
{
  static constexpr int NFNC = PairBodyRoundedPolyhedron::NFNC;
  Atom *atom = this->atom;
  Memory *memory = this->memory;

  this->ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = this->comm->nthreads;
  const int inum = this->list->inum;

  // grow the per-atom lists if necessary and initialize

  if (atom->nmax > this->nmax) {
    memory->destroy(this->dnum);
    memory->destroy(this->dfirst);
    memory->destroy(this->ednum);
    memory->destroy(this->edfirst);
    memory->destroy(this->facnum);
    memory->destroy(this->facfirst);
    memory->destroy(this->enclosing_radius);
    memory->destroy(this->rounded_radius);
    this->nmax = atom->nmax;
    memory->create(this->dnum, this->nmax, "pair:dnum");
    memory->create(this->dfirst, this->nmax, "pair:dfirst");
    memory->create(this->ednum, this->nmax, "pair:ednum");
    memory->create(this->edfirst, this->nmax, "pair:edfirst");
    memory->create(this->facnum, this->nmax, "pair:facnum");
    memory->create(this->facfirst, this->nmax, "pair:facfirst");
    memory->create(this->enclosing_radius, this->nmax, "pair:enclosing_radius");
    memory->create(this->rounded_radius, this->nmax, "pair:rounded_radius");
  }

  this->ndiscrete = this->nedge = this->nface = 0;
  for (int i = 0; i < nall; i++) this->dnum[i] = this->ednum[i] = this->facnum[i] = 0;

  // body2space modifies the shared vertex, edge and face lists,
  // so it is called for all bodies before the parallel region

  int *body = atom->body;
  for (int i = 0; i < nall; i++)
    if (body[i] >= 0) this->body2space(i);

  // per-thread forces and torques that do not derive from the energy,
  // per-thread scratch space, and per-thread data of the contact models

  if (atom->nmax > this->nmax_fnc) {
    memory->destroy(this->fnc);
    this->nmax_fnc = atom->nmax;
    memory->create(this->fnc, this->nmax_fnc, NFNC, "pair:fnc");
  }
  if ((atom->nmax > nmax_thr) || ((int) scratch_thr.size() != nthreads)) {
    memory->destroy(fnc_thr);
    nmax_thr = atom->nmax;
    memory->create(fnc_thr, nthreads * nmax_thr, NFNC, "pair:fnc_thr");
    scratch_thr.resize(nthreads);
  }
  init_thr(nthreads);

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, this->eatom, this->vatom, nullptr, thr);

    double **fnc_t = fnc_thr + tid * nmax_thr;
    for (int i = 0; i < nall; i++)
      for (int k = 0; k < NFNC; k++) fnc_t[i][k] = 0.0;

    if (this->evflag) {
      if (eflag) {
        if (this->force->newton_pair) eval<1, 1, 1>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
        else eval<1, 1, 0>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
      } else {
        if (this->force->newton_pair) eval<1, 0, 1>(ifrom, ito, thr, fnc_t, scratch_thr[tid]);
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
      this->fnc[i][k] = sum;
    }
  }

  this->work_nonconservative();
}

/* ---------------------------------------------------------------------- */

template <class BASE>
template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBodyRoundedPolyhedronOMPT<BASE>::eval(int iifrom, int iito, ThrData *const thr,
                                                double **fnc_t,
                                                PairBodyRoundedPolyhedron::Scratch &s)
{
  Atom *atom = this->atom;
  double **x = atom->x;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double **f = thr->get_f();
  double **torque = thr->get_torque();
  double *radius = atom->radius;
  int *body = atom->body;
  const int nlocal = atom->nlocal;
  const int history = this->history;
  const int size_history = this->size_history;
  const double cut_inner = this->cut_inner;

  const int *const ilist = this->list->ilist;
  const int *const numneigh = this->list->numneigh;
  int **const firstneigh = this->list->firstneigh;

  double facc[3];

  // tangential displacements of the pairs in the neighbor list,
  // each thread updates those of its own atoms i

  int *touch = nullptr;
  double *allshear = nullptr;
  s.shear = nullptr;

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const double radi = radius[i];
    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    if (history) {
      touch = this->fix_history->firstflag[i];
      allshear = this->fix_history->firstvalue[i];
    }

    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj] & NEIGHMASK;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;

      // the tangential displacement is reset unless the pair is in contact

      if (history) {
        s.shear = &allshear[size_history * jj];
        s.shear_i = i;
        s.touched = 0;
        s.fnsum = 0.0;
        for (int k = 0; k < 3; k++) s.pcsum[k] = s.vtsum[k] = s.nsum[k] = 0.0;
      }

      // no interaction

      if ((body[i] >= 0) && (body[j] >= 0) && (sqrt(rsq) <= radi + radius[j] + cut_inner)) {
        double evdwl = 0.0;
        facc[0] = facc[1] = facc[2] = 0.0;
        this->pair_interaction(i, j, delx, dely, delz, rsq, x, v, angmom, f, torque, fnc_t, s,
                               evdwl, facc);
        if (history) this->history_friction(i, j, x, f, torque, fnc_t, s, facc);

        if (EVFLAG)
          ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR, EFLAG ? evdwl : 0.0, 0.0, facc[0],
                           facc[1], facc[2], delx, dely, delz, thr);
      }

      if (history) {
        touch[jj] = s.touched;
        if (!s.touched)
          for (int k = 0; k < size_history; k++) s.shear[k] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class BASE> double PairBodyRoundedPolyhedronOMPT<BASE>::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += BASE::memory_usage();
  bytes += (double) nmax_thr * scratch_thr.size() * PairBodyRoundedPolyhedron::NFNC *
      sizeof(double);
  return bytes;
}

}    // namespace LAMMPS_NS

#endif
