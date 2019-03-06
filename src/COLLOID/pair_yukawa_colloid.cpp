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
   Contributing authors: Randy Schunk (Sandia)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include "pair_yukawa_colloid.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairYukawaColloid::PairYukawaColloid(LAMMPS *lmp) : PairYukawa(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

void PairYukawaColloid::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,radi,radj;
  double rsq,r,rinv,screening,forceyukawa,factor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        screening = exp(-kappa*(r-(radi+radj)));
        forceyukawa = a[itype][jtype] * screening;

        fpair = factor*forceyukawa * rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = a[itype][jtype]/kappa * screening - offset[itype][jtype];
          evdwl *= factor;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairYukawaColloid::init_style()
{
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair yukawa/colloid requires atom style sphere");

  neighbor->request(this,instance_me);

  // require that atom radii are identical within each type

  for (int i = 1; i <= atom->ntypes; i++)
    if (!atom->radius_consistency(i,rad[i]))
      error->all(FLERR,"Pair yukawa/colloid requires atoms with same type "
                 "have same radius");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairYukawaColloid::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    a[i][j] = mix_energy(a[i][i],a[j][j],1.0,1.0);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (offset_flag && (kappa != 0.0)) {
    double screening = exp(-kappa * (cut[i][j] - (rad[i]+rad[j])));
    offset[i][j] = a[i][j]/kappa * screening;
  } else offset[i][j] = 0.0;

  a[j][i] = a[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairYukawaColloid::single(int /*i*/, int /*j*/, int itype, int jtype,
                                 double rsq,
                                 double /*factor_coul*/, double factor_lj,
                                 double &fforce)
{
  double r,rinv,screening,forceyukawa,phi;

  r = sqrt(rsq);
  rinv = 1.0/r;
  screening = exp(-kappa*(r-(rad[itype]+rad[jtype])));
  forceyukawa = a[itype][jtype] * screening;
  fforce = factor_lj*forceyukawa * rinv;

  phi = a[itype][jtype]/kappa * screening  - offset[itype][jtype];
  return factor_lj*phi;
}
