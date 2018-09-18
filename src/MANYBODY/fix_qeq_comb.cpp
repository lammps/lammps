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
   Contributing authors: Ray Shan (Sandia, tnshan@sandia.gov)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pair_comb.h"
#include "pair_comb3.h"
#include "fix_qeq_comb.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "group.h"
#include "respa.h"
#include "update.h"
#include "memory.h"
#include "error.h"
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixQEQComb::FixQEQComb(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
  fp(NULL), comb(NULL), comb3(NULL), qf(NULL), q1(NULL), q2(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal fix qeq/comb command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  nevery = force->inumeric(FLERR,arg[3]);
  precision = force->numeric(FLERR,arg[4]);

  if (nevery <= 0 || precision <= 0.0)
    error->all(FLERR,"Illegal fix qeq/comb command");

  MPI_Comm_rank(world,&me);

  // optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/comb command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          snprintf(str,128,"Cannot open fix qeq/comb file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qeq/comb command");
  }

  nmax = atom->nmax;
  memory->create(qf,nmax,"qeq:qf");
  memory->create(q1,nmax,"qeq:q1");
  memory->create(q2,nmax,"qeq:q2");
  vector_atom = qf;

  // zero the vector since dump may access it on timestep 0
  // zero the vector since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) qf[i] = 0.0;

  comb = NULL;
  comb3 = NULL;

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixQEQComb::~FixQEQComb()
{
  if (me == 0 && fp) fclose(fp);
  memory->destroy(qf);
  memory->destroy(q1);
  memory->destroy(q2);
}

/* ---------------------------------------------------------------------- */

int FixQEQComb::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQEQComb::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/comb requires atom attribute q");

  comb = (PairComb *) force->pair_match("comb",1);
  comb3 = (PairComb3 *) force->pair_match("comb3",1);
  if (comb == NULL && comb3 == NULL)
    error->all(FLERR,"Must use pair_style comb or comb3 with fix qeq/comb");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/comb group has no atoms");
}

/* ---------------------------------------------------------------------- */

void FixQEQComb::setup(int vflag)
{
  firstflag = 1;
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
  firstflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixQEQComb::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEQComb::post_force(int /*vflag*/)
{
  int i,ii,iloop,loopmax,inum,*ilist;
  double heatpq,qmass,dtq,dtq2;
  double enegchkall,enegmaxall;

  if (update->ntimestep % nevery) return;

  // reallocate work arrays if necessary
  // qf = charge force
  // q1 = charge displacement
  // q2 = tmp storage of charge force for next iteration

  if (atom->nmax > nmax) {
    memory->destroy(qf);
    memory->destroy(q1);
    memory->destroy(q2);
    nmax = atom->nmax;
    memory->create(qf,nmax,"qeq:qf");
    memory->create(q1,nmax,"qeq:q1");
    memory->create(q2,nmax,"qeq:q2");
    vector_atom = qf;
  }

  // more loops for first-time charge equilibrium

  iloop = 0;
  if (firstflag) loopmax = 200;
  else loopmax = 100;

  // charge-equilibration loop

  if (me == 0 && fp)
    fprintf(fp,"Charge equilibration on step " BIGINT_FORMAT "\n",
            update->ntimestep);

  heatpq = 0.05;
  qmass  = 0.016;
  dtq    = 0.01;
  dtq2   = 0.5*dtq*dtq/qmass;

  double enegchk = 0.0;
  double enegtot = 0.0;
  double enegmax = 0.0;

  double *q = atom->q;
  int *mask = atom->mask;

 if (comb) {
    inum = comb->list->inum;
    ilist = comb->list->ilist;
  }
  if (comb3) {
    inum = comb3->list->inum;
    ilist = comb3->list->ilist;
  }
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    q1[i] = q2[i] = qf[i] = 0.0;
  }

  for (iloop = 0; iloop < loopmax; iloop ++ ) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        q1[i] += qf[i]*dtq2 - heatpq*q1[i];
        q[i]  += q1[i];
      }
    }

    comm->forward_comm_fix(this);
    if(comb) enegtot = comb->yasu_char(qf,igroup);
    if(comb3) enegtot = comb3->combqeq(qf,igroup);

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

    if (enegchk <= precision && enegmax <= 100.0*precision) break;

    if (me == 0 && fp)
      fprintf(fp,"  iteration: %d, enegtot %.6g, "
              "enegmax %.6g, fq deviation: %.6g\n",
              iloop,enegtot,enegmax,enegchk);

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit)
        q1[i] += qf[i]*dtq2 - heatpq*q1[i];
    }
  }

  if (me == 0 && fp) {
    if (iloop == loopmax)
      fprintf(fp,"Charges did not converge in %d iterations\n",iloop);
    else
      fprintf(fp,"Charges converged in %d iterations to %.10f tolerance\n",
              iloop,enegchk);
  }
}

/* ---------------------------------------------------------------------- */

void FixQEQComb::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEQComb::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}
/* ---------------------------------------------------------------------- */

int FixQEQComb::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i ++) {
    j = list[i];
    buf[m++] = atom->q[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEQComb::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  for (i = first; i < last; i++) atom->q[i] = buf[m++];
}
