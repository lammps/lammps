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
   Contributing authors:
     James Fischer, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natoli, Stone Ridge Technology
------------------------------------------------------------------------- */

#include "pair_morse_opt.h"
#include <cmath>

#include "atom.h"
#include "force.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMorseOpt::PairMorseOpt(LAMMPS *lmp) : PairMorse(lmp) {}

/* ---------------------------------------------------------------------- */

void PairMorseOpt::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) return eval<1,1,1>();
      else return eval<1,1,0>();
    } else {
      if (force->newton_pair) return eval<1,0,1>();
      else return eval<1,0,0>();
    }
  } else {
    if (force->newton_pair) return eval<0,0,1>();
    else return eval<0,0,0>();
  }
}

/* ---------------------------------------------------------------------- */

template < int EVFLAG, int EFLAG, int NEWTON_PAIR >
void PairMorseOpt::eval()
{
  typedef struct { double x,y,z; } vec3_t;

  typedef struct {
    double cutsq,r0,alpha,morse1,d0,offset;
    double _pad[2];
  } fast_alpha_t;

  int i,j,ii,jj,inum,jnum,itype,jtype,sbindex;
  double factor_lj;
  double evdwl = 0.0;

  double** _noalias x = atom->x;
  double** _noalias f = atom->f;
  int* _noalias type = atom->type;
  int nlocal = atom->nlocal;
  double* _noalias special_lj = force->special_lj;

  inum = list->inum;
  int* _noalias ilist = list->ilist;
  int** _noalias firstneigh = list->firstneigh;
  int* _noalias numneigh = list->numneigh;

  auto * _noalias xx = (vec3_t*)x[0];
  auto * _noalias ff = (vec3_t*)f[0];

  int ntypes = atom->ntypes;
  int ntypes2 = ntypes*ntypes;

  auto * _noalias fast_alpha =
    (fast_alpha_t*) malloc(ntypes2*sizeof(fast_alpha_t));
  for (i = 0; i < ntypes; i++) for (j = 0; j < ntypes; j++) {
    fast_alpha_t& a = fast_alpha[i*ntypes+j];
    a.cutsq = cutsq[i+1][j+1];
    a.r0 = r0[i+1][j+1];
    a.alpha = alpha[i+1][j+1];
    a.morse1 = morse1[i+1][j+1];
    a.d0 = d0[i+1][j+1];
    a.offset = offset[i+1][j+1];
  }
  auto * _noalias tabsix = fast_alpha;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    itype = type[i] - 1;
    int* _noalias jlist = firstneigh[i];
    jnum = numneigh[i];

    double tmpfx = 0.0;
    double tmpfy = 0.0;
    double tmpfz = 0.0;

    auto * _noalias tabsixi = (fast_alpha_t*)&tabsix[itype*ntypes];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      sbindex = sbmask(j);

      if (sbindex == 0) {
        double delx = xtmp - xx[j].x;
        double dely = ytmp - xx[j].y;
        double delz = ztmp - xx[j].z;
        double rsq = delx*delx + dely*dely + delz*delz;

        jtype = type[j] - 1;

        fast_alpha_t& a = tabsixi[jtype];
        if (rsq < a.cutsq) {
          double r = sqrt(rsq);
          double dr = r - a.r0;
          double dexp = exp(-a.alpha * dr);
          double fpair = a.morse1 * (dexp*dexp - dexp) / r;

          tmpfx += delx*fpair;
          tmpfy += dely*fpair;
          tmpfz += delz*fpair;
          if (NEWTON_PAIR || j < nlocal) {
            ff[j].x -= delx*fpair;
            ff[j].y -= dely*fpair;
            ff[j].z -= delz*fpair;
          }

          if (EFLAG) evdwl = a.d0 * (dexp*dexp - 2.0*dexp) - a.offset;

          if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                               evdwl,0.0,fpair,delx,dely,delz);
        }

      } else {
        factor_lj = special_lj[sbindex];
        j &= NEIGHMASK;

        double delx = xtmp - xx[j].x;
        double dely = ytmp - xx[j].y;
        double delz = ztmp - xx[j].z;
        double rsq = delx*delx + dely*dely + delz*delz;

        jtype = type[j] - 1;

        fast_alpha_t& a = tabsixi[jtype];
        if (rsq < a.cutsq) {
          double r = sqrt(rsq);
          double dr = r - a.r0;
          double dexp = exp(-a.alpha * dr);
          double fpair = factor_lj * a.morse1 * (dexp*dexp - dexp) / r;

          tmpfx += delx*fpair;
          tmpfy += dely*fpair;
          tmpfz += delz*fpair;
          if (NEWTON_PAIR || j < nlocal) {
            ff[j].x -= delx*fpair;
            ff[j].y -= dely*fpair;
            ff[j].z -= delz*fpair;
          }

          if (EFLAG) {
            evdwl = a.d0 * (dexp*dexp - 2.0*dexp) - a.offset;
            evdwl *= factor_lj;
          }

          if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                               evdwl,0.0,fpair,delx,dely,delz);
        }
      }
    }

    ff[i].x += tmpfx;
    ff[i].y += tmpfy;
    ff[i].z += tmpfz;
  }

  free(fast_alpha); fast_alpha = nullptr;

  if (vflag_fdotr) virial_fdotr_compute();
}
