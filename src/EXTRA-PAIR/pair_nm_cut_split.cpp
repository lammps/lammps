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
   Contributing Author: Julien Devemy (ICCF), Robert S. Hoy (USF), Joseph D. Dietz (USF)
------------------------------------------------------------------------- */

#include "pair_nm_cut_split.h"

#include "atom.h"
#include "force.h"
#include "math_special.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathSpecial::powint;

/* ---------------------------------------------------------------------- */
PairNMCutSplit::PairNMCutSplit(LAMMPS *lmp) : PairNMCut(lmp)
{
  writedata = 1;
}

void PairNMCutSplit::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,factor_lj;
  double r,forcenm,rminv,rninv;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
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
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);

        // r < r0 --> use generalized LJ
        if (rsq < r0[itype][jtype]*r0[itype][jtype]) {
          forcenm = e0nm[itype][jtype]*nm[itype][jtype]*
            (r0n[itype][jtype]/pow(r,nn[itype][jtype])
             -r0m[itype][jtype]/pow(r,mm[itype][jtype]));
        }
        // r > r0 --> use standard LJ (m = 6 n = 12)
        else forcenm =(e0[itype][jtype]/6.0)*72.0*(4.0/powint(r,12)-2.0/powint(r,6));

        fpair = factor_lj*forcenm*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          // r < r0 --> use generalized LJ
          if (rsq < r0[itype][jtype]*r0[itype][jtype]) {
            rminv = pow(r2inv,mm[itype][jtype]/2.0);
            rninv = pow(r2inv,nn[itype][jtype]/2.0);

            evdwl = e0nm[itype][jtype]*(mm[itype][jtype]*r0n[itype][jtype]*rninv -
                                        nn[itype][jtype]*r0m[itype][jtype]*rminv) -
              offset[itype][jtype];
          }
          // r > r0 --> use standard LJ (m = 6 n = 12)
          else evdwl = (e0[itype][jtype]/6.0)*(24.0*powint(r2inv,6) - 24.0*pow(r2inv,3));
          evdwl *= factor_lj;
        }
        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairNMCutSplit::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq, double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r2inv,rminv,rninv,r,forcenm,phinm;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rminv = pow(r2inv,mm[itype][jtype]/2.0);
  rninv = pow(r2inv,nn[itype][jtype]/2.0);
  // r < r_0, use generalized LJ
  if (rsq < r0[itype][jtype]*r0[itype][jtype]) {  // note the addition of the r0 factor
     forcenm = e0nm[itype][jtype]*nm[itype][jtype]*
      (r0n[itype][jtype]/pow(r,nn[itype][jtype])-r0m[itype][jtype]/pow(r,mm[itype][jtype]));
      phinm = e0nm[itype][jtype]*(mm[itype][jtype]*r0n[itype][jtype]*rninv
      -nn[itype][jtype]*r0m[itype][jtype]*rminv)-offset[itype][jtype];

  }
  // r > r_0 --> use standard LJ (m = 6 n = 12)
  else {
    forcenm = (e0[itype][jtype]/6.0)*72.0*(4.0/powint(r,12)-2.0/powint(r,6));
    phinm = (e0[itype][jtype]/6.0)*(24.0*powint(r2inv,6)-24.0*powint(r2inv,3));
  }

  fforce = factor_lj*forcenm*r2inv;
  return factor_lj*phinm;
}

