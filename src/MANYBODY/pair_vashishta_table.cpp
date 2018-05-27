/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Anders Hafreager (UiO), andershaf@gmail.com
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_vashishta_table.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairVashishtaTable::PairVashishtaTable(LAMMPS *lmp) : PairVashishta(lmp)
{
  forceTable = NULL;
  potentialTable = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairVashishtaTable::~PairVashishtaTable()
{
  memory->destroy(forceTable);
  memory->destroy(potentialTable);
}

/* ---------------------------------------------------------------------- */

void PairVashishtaTable::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const double cutshortsq = r0max*r0max;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutshortsq) {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      if (rsq >= params[ijparam].cutsq) continue;

      twobody_table(params[ijparam],rsq,fpair,eflag,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                   evdwl,0.0,fpair,delx,dely,delz);
    }

    jnumm1 = numshort - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[ijparam].cutsq2) continue;

      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        ktype = map[type[k]];
        ikparam = elem2param[itype][ktype][ktype];
        ijkparam = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[ikparam].cutsq2) continue;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

        fxtmp -= fj[0] + fk[0];
        fytmp -= fj[1] + fk[1];
        fztmp -= fj[2] + fk[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairVashishtaTable::twobody_table(const Param &param, double rsq,
                                       double &fforce, int eflag, double &eng)
{
  // use analytic form if rsq is inside inner cutoff

  if (rsq < tabinnersq) {
    Param *pparam = const_cast<Param *> (&param);
    PairVashishta::twobody(pparam,rsq,fforce,eflag,eng);
    return;
  }

  // double -> int will only keep the 0.xxxx part

  const int tableIndex = (rsq - tabinnersq)*oneOverDeltaR2;
  const double fraction = (rsq - tabinnersq)*oneOverDeltaR2 - tableIndex;

  // force/energy are linearly interpolated between two adjacent values

  double force0 = forceTable[param.ielement][param.jelement][tableIndex];
  double force1 = forceTable[param.ielement][param.jelement][tableIndex+1];
  fforce = (1.0 - fraction)*force0 + fraction*force1;

  if (evflag) {
    double energy0 = potentialTable[param.ielement][param.jelement][tableIndex];
    double energy1 = potentialTable[param.ielement][param.jelement]
      [tableIndex+1];
    eng = (1.0 - fraction)*energy0 + fraction*energy1;
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairVashishtaTable::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  ntable = force->inumeric(FLERR,arg[0]);
  tabinner = force->numeric(FLERR,arg[1]);

  if (tabinner <= 0.0)
    error->all(FLERR,"Illegal inner cutoff for tabulation");
}

/* ---------------------------------------------------------------------- */

void PairVashishtaTable::setup_params()
{
  PairVashishta::setup_params();

  create_tables();
}

/* ---------------------------------------------------------------------- */

void PairVashishtaTable::create_tables()
{
  memory->destroy(forceTable);
  memory->destroy(potentialTable);
  forceTable = NULL;
  potentialTable = NULL;

  tabinnersq = tabinner*tabinner;

  deltaR2 = (cutmax*cutmax - tabinnersq) / (ntable-1);
  oneOverDeltaR2 = 1.0/deltaR2;

  memory->create(forceTable,nelements,nelements,ntable+1,
                 "pair:vashishta:forceTable");
  memory->create(potentialTable,nelements,nelements,ntable+1,
                 "pair:vashishta:potentialTable");

  // tabulalate energy/force via analytic twobody() in parent

  int i,j,idx;
  double rsq,fpair,eng;

  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {
      int ijparam = elem2param[i][j][j];
      for (idx = 0; idx <= ntable; idx++) {
        rsq = tabinnersq + idx*deltaR2;
        PairVashishta::twobody(&params[ijparam],rsq,fpair,1,eng);
        forceTable[i][j][idx] = fpair;
        potentialTable[i][j][idx] = eng;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of tabulation arrays
------------------------------------------------------------------------- */

double PairVashishtaTable::memory_usage()
{
  double bytes = 2*nelements*nelements*sizeof(double)*ntable;
  return bytes;
}
