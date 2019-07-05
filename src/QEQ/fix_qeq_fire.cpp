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
   Contributing author: Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "fix_qeq_fire.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "pair.h"
#include "pair_comb.h"
#include "pair_comb3.h"
#include "kspace.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELAYSTEP 0
#define DT_GROW 1.1
#define DT_SHRINK 0.5
#define ALPHA0 0.8
#define ALPHA_SHRINK 0.10
#define TMAX 10.0

/* ---------------------------------------------------------------------- */

FixQEqFire::FixQEqFire(LAMMPS *lmp, int narg, char **arg) :
  FixQEq(lmp, narg, arg)
{
  qdamp = 0.20;
  qstep = 0.20;

  int iarg = 8;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"qdamp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/fire command");
      qdamp = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"qstep") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/fire command");
      qstep = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qeq/fire command");
  }

  comb = NULL;
  comb3 = NULL;
}

/* ---------------------------------------------------------------------- */

void FixQEqFire::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/fire requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/fire group has no atoms");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;

  if (tolerance < 1e-4)
    if (comm->me == 0)
      error->warning(FLERR,"Fix qeq/fire tolerance may be too small"
                    " for damped fires");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  comb = (PairComb *) force->pair_match("comb",1);
  comb3 = (PairComb3 *) force->pair_match("comb3",1);

}

/* ---------------------------------------------------------------------- */

void FixQEqFire::pre_force(int /*vflag*/)
{
  int inum, *ilist;
  int i,ii,iloop;

  double *q = atom->q;
  double vmax,vdotf,vdotfall,vdotv,vdotvall,fdotf,fdotfall;
  double scale1,scale2;
  double dtvone,dtv;
  double enegtot,enegchk;
  double alpha = qdamp;
  double dt, dtmax;
  double enegchkall;
  bigint ntimestep = update->ntimestep;
  bigint last_negative = 0;

  if (ntimestep % nevery) return;

  if (atom->nmax > nmax) reallocate_storage();

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qv[i] = 0.0;
  }

  dt = qstep;
  dtmax = TMAX * dt;

  for (iloop = 0; iloop < maxiter; iloop ++ ) {
    pack_flag = 1;
    comm->forward_comm_fix(this);

    if (comb) {
      comb->yasu_char(qf,igroup);
      enegtot = comb->enegtot / ngroup;
    } else if (comb3) {
      comb3->combqeq(qf,igroup);
      enegtot = comb3->enegtot / ngroup;
    } else {
      enegtot = compute_eneg();
      enegtot /= ngroup;
    }

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      qf[i] -= enegtot;         // Enforce adiabatic
    }

    // FIRE minimization algorithm
    // vdotfall = v dot f = qv dot qf
    vdotf = 0.0;
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      vdotf += (qv[i]*qf[i]);
    }
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    if (vdotfall > 0.0) {
      vdotv = fdotf = 0.0;
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        vdotv += qv[i]*qv[i];
        fdotf += qf[i]*qf[i];
      }
      MPI_Allreduce(&vdotv,&vdotvall,1,MPI_DOUBLE,MPI_SUM,world);
      MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,world);

      scale1 = 1.0 - alpha;
      if (fdotfall == 0.0) scale2 = 0.0;
      else scale2 = alpha * sqrt(vdotvall/fdotfall);

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        qv[i] = scale1*qv[i] + scale2*qf[i];
      }
      if (ntimestep - last_negative > DELAYSTEP) {
        dt = MIN(dt*DT_GROW,dtmax);
        alpha *= ALPHA_SHRINK;
      }
    } else {
      last_negative = ntimestep;
      dt *= DT_SHRINK;
      alpha = ALPHA0;
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        qv[i] = 0.0;
      }
    }

    // limit timestep so no charges change more than dmax
    dtvone = dt;
    double dmax = 0.1;
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      vmax = MAX(fabs(qv[i]),0);
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }
    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);
    //dtv = dt;

    // Euler integration step
    enegchk = 0.0;
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      q[i] -= dtv * qv[i];
      qv[i] += dtv * qf[i];
      enegchk += fabs(qf[i]);
    }
    MPI_Allreduce(&enegchk,&enegchkall,1,MPI_DOUBLE,MPI_SUM,world);
    enegchk = enegchkall / ngroup;

    if (enegchk < tolerance) break;
  }

  if (comm->me == 0) {
    if (iloop == maxiter) {
      char str[128];
      sprintf(str,"Charges did not converge at step " BIGINT_FORMAT
                  ": %lg",update->ntimestep,enegchk);
      error->warning(FLERR,str);
    }
  }

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

double FixQEqFire::compute_eneg()
{
  int i, j, ii, jj, inum, jnum, itype;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double eneg, enegtot;
  double r, rsq, delr[3], rinv;

  int *type = atom->type;
  int *mask = atom->mask;
  double *q = atom->q;
  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      qf[i] = 0.0;
  }

  // communicating charge force to all nodes, first forward then reverse
  pack_flag = 2;
  comm->forward_comm_fix(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    if (mask[i] & groupbit) {

      qf[i] += chi[itype] + eta[itype] * q[i];

      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delr[0] = x[i][0] - x[j][0];
        delr[1] = x[i][1] - x[j][1];
        delr[2] = x[i][2] - x[j][2];
        rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];

        if (rsq > cutoff_sq) continue;

        r = sqrt(rsq);
        rinv = 1.0/r;
        qf[i] += q[j] * rinv;
        qf[j] += q[i] * rinv;
      }
    }
  }

  pack_flag = 2;
  comm->reverse_comm_fix(this);

  // sum charge force on each node and return it

  eneg = enegtot = 0.0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      eneg += qf[i];
  }
  MPI_Allreduce(&eneg,&enegtot,1,MPI_DOUBLE,MPI_SUM,world);
  return enegtot;

}

/* ---------------------------------------------------------------------- */

int FixQEqFire::pack_forward_comm(int n, int *list, double *buf,
                          int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;

  if (pack_flag == 1)
    for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  else if( pack_flag == 2 )
    for(m = 0; m < n; m++) buf[m] = qf[list[m]];

  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEqFire::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
  else if( pack_flag == 2)
    for(m = 0, i = first; m < n; m++, i++) qf[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

int FixQEqFire::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  for(m = 0, i = first; m < n; m++, i++) buf[m] = qf[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEqFire::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m;

  for(m = 0; m < n; m++) qf[list[m]] += buf[m];
}

/* ---------------------------------------------------------------------- */
