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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_sw_table.h"
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

PairSWTable::PairSWTable(LAMMPS *lmp) : PairSW(lmp)
{
  neigh3BodyMax = 0;
  neigh3BodyCount = NULL; 
  neigh3Body = NULL;
  forceTable = NULL;
  potentialTable = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairSWTable::~PairSWTable()
{
  memory->destroy(forceTable);
  memory->destroy(potentialTable);
  memory->destroy(neigh3BodyCount);
  memory->destroy(neigh3Body);
}

/* ---------------------------------------------------------------------- */

void PairSWTable::compute(int eflag, int vflag)
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

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // reallocate 3-body neighbor list if necessary
  // NOTE: using 1000 is inefficient
  //       could make this a LAMMPS paged neighbor list

  if (nlocal > neigh3BodyMax) {
    neigh3BodyMax = atom->nmax;
    memory->destroy(neigh3BodyCount);
    memory->destroy(neigh3Body);
    memory->create(neigh3BodyCount,neigh3BodyMax,
                   "pair:sw:neigh3BodyCount");
    memory->create(neigh3Body,neigh3BodyMax,1000,
                   "pair:sw:neigh3Body");
  }

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // reset the 3-body neighbor list

    neigh3BodyCount[i] = 0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];

      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      ijparam = elem2param[itype][jtype][jtype];

      if (rsq < params[ijparam].cutsq) {
        neigh3Body[i][neigh3BodyCount[i]] = j;
        neigh3BodyCount[i]++;
      }

      if (rsq >= params[ijparam].cutsq) continue;

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }      
      twobody_table(params[ijparam],rsq,fpair,eflag,evdwl);

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    jlist = neigh3Body[i];
    jnum = neigh3BodyCount[i];
    jnumm1 = jnum - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[ijparam].cutsq) continue;

      for (kk = jj+1; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem2param[itype][ktype][ktype];
        ijkparam = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[ikparam].cutsq) continue;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

        f[i][0] -= fj[0] + fk[0];
        f[i][1] -= fj[1] + fk[1];
        f[i][2] -= fj[2] + fk[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */
void PairSWTable::twobody_table(Param &param, double rsq, 
                                       double &fforce, int eflag, double &eng)
{
  // use analytic form if rsq is inside inner cutoff

  if (rsq < tabinnersq) {
    Param *pparam = const_cast<Param *> (&param);
    PairSW::twobody(pparam,rsq,fforce,eflag,eng);
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
    double energy1 = potentialTable[param.ielement][param.jelement][tableIndex+1];
    eng = (1.0 - fraction)*energy0 + fraction*energy1;
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSWTable::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  ntable = force->inumeric(FLERR,arg[0]);
  tabinner = force->numeric(FLERR,arg[1]);

  if (tabinner <= 0.0)
    error->all(FLERR,"Illegal inner cutoff for tabulation");
}

/* ---------------------------------------------------------------------- */

void PairSWTable::setup_params()
{
  PairSW::setup_params();

  create_tables();
}

/* ---------------------------------------------------------------------- */

void PairSWTable::create_tables()
{
  memory->destroy(forceTable);
  memory->destroy(potentialTable);
  forceTable = NULL;
  potentialTable = NULL;

  tabinnersq = tabinner*tabinner;

  deltaR2 = (cutmax*cutmax - tabinnersq) / (ntable-1);
  oneOverDeltaR2 = 1.0/deltaR2;

  memory->create(forceTable,nelements,nelements,ntable+1,
                 "pair:sw:forceTable");
  memory->create(potentialTable,nelements,nelements,ntable+1,
                 "pair:sw:potentialTable");

  // tabulalate energy/force via analytic twobody() in parent

  int i,j,idx;
  double rsq,fpair,eng;

  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {
      int ijparam = elem2param[i][j][j];
      for (idx = 0; idx <= ntable; idx++) {
        rsq = tabinnersq + idx*deltaR2;
        PairSW::twobody(&params[ijparam],rsq,fpair,1,eng);
        if(rsq >= params[ijparam].cutsq) {
          // the rainv factor might blow up for rsq > cutsq but that isn't physical. Force and energy goes to zero here.
          fpair = 0.0;
          eng = 0.0;
        }

        forceTable[i][j][idx] = fpair;
        potentialTable[i][j][idx] = eng;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of tabulation arrays
------------------------------------------------------------------------- */

double PairSWTable::memory_usage()
{
  double bytes = 2*nelements*nelements*sizeof(double)*ntable;
  return bytes;
}
