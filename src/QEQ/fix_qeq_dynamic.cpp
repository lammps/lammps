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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_qeq_dynamic.h"
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
#include "kspace.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqDynamic::FixQEqDynamic(LAMMPS *lmp, int narg, char **arg) : 
  FixQEq(lmp, narg, arg) 
{
  qdamp = 0.10;
  qstep = 0.02;

  int iarg = 8;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"qdamp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/dynamic command");
      qdamp = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"qstep") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/dynamic command");
      qstep = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qeq/dynamic command");
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqDynamic::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/dynamic requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/dynamic group has no atoms");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;

  if (tolerance < 1e-4) 
    if (comm->me == 0)
      error->warning(FLERR,"Fix qeq/dynamic tolerance may be too small"
		    " for damped dynamics");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixQEqDynamic::pre_force(int vflag)
{
  int i,ii,iloop,inum,*ilist;
  double qmass,dtq2;
  double enegchkall,enegmaxall;

  double *q = atom->q;
  int *mask = atom->mask;

  double enegchk = 0.0;
  double enegtot = 0.0;
  double enegmax = 0.0;

  if (update->ntimestep % nevery) return;

  n = atom->nlocal;
  N = atom->nlocal + atom->nghost;

  if( atom->nmax > nmax ) reallocate_storage();

  inum = list->inum;
  ilist = list->ilist;

  qmass  = 0.016;
  dtq2   = 0.5*qstep*qstep/qmass;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    q1[i] = q2[i] = qf[i] = 0.0;
  }

  for (iloop = 0; iloop < maxiter; iloop ++ ) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        q1[i] += qf[i]*dtq2 - qdamp*q1[i];
        q[i]  += q1[i];
      }
    }

    pack_flag = 1;
    comm->forward_comm_fix(this);

    enegtot = compute_eneg();
    enegtot /= ngroup;

    enegchk = enegmax = 0.0;

    for (ii = 0; ii < inum ; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        q2[i] = enegtot-qf[i];
        enegmax = MAX(enegmax,fabs(q2[i]));
        enegchk += fabs(q2[i]);
        qf[i] = q2[i];
      }
    }

    MPI_Allreduce(&enegchk,&enegchkall,1,MPI_DOUBLE,MPI_SUM,world);
    enegchk = enegchkall/ngroup;
    MPI_Allreduce(&enegmax,&enegmaxall,1,MPI_DOUBLE,MPI_MAX,world);
    enegmax = enegmaxall;

    if (enegchk <= tolerance && enegmax <= 100.0*tolerance) break;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit)
        q1[i] += qf[i]*dtq2 - qdamp*q1[i];
    }
  }

  if (comm->me == 0) {
    if (iloop == maxiter) {
      char str[128];
      sprintf(str,"Charges did not converge at step "BIGINT_FORMAT
		  ": %lg",update->ntimestep,enegchk);
      error->warning(FLERR,str);
    }
  }

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

double FixQEqDynamic::compute_eneg()
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

int FixQEqDynamic::pack_forward_comm(int n, int *list, double *buf,
                          int pbc_flag, int *pbc)
{
  int m;

  if( pack_flag == 1 )
    for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  else if( pack_flag == 2 )
    for(m = 0; m < n; m++) buf[m] = qf[list[m]];

  return n;
}

/* ---------------------------------------------------------------------- */

void FixQEqDynamic::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if( pack_flag == 1)
    for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
  else if( pack_flag == 2)
    for(m = 0, i = first; m < n; m++, i++) qf[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

int FixQEqDynamic::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  for(m = 0, i = first; m < n; m++, i++) buf[m] = qf[i];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixQEqDynamic::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m;

  for(m = 0; m < n; m++) qf[list[m]] += buf[m];
}

/* ---------------------------------------------------------------------- */
