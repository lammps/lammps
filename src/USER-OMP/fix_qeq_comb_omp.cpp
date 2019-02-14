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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstring>
#include "fix_qeq_comb_omp.h"
#include "fix_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "respa.h"
#include "update.h"
#include "pair_comb_omp.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixQEQCombOMP::FixQEQCombOMP(LAMMPS *lmp, int narg, char **arg) :
  FixQEQComb(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix qeq/comb/omp command");
}

/* ---------------------------------------------------------------------- */

void FixQEQCombOMP::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/comb/omp requires atom attribute q");

  if (NULL != force->pair_match("comb3",0))
    error->all(FLERR,"No support for comb3 currently available in USER-OMP");

  comb = (PairComb *) force->pair_match("comb/omp",1);
  if (comb == NULL)
    comb = (PairComb *) force->pair_match("comb",1);
  if (comb == NULL)
    error->all(FLERR,"Must use pair_style comb or "
               "comb/omp with fix qeq/comb/omp");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/comb group has no atoms");
}

/* ---------------------------------------------------------------------- */

void FixQEQCombOMP::post_force(int /* vflag */)
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
  if (firstflag) loopmax = 500;
  else loopmax = 200;

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

  inum = comb->list->inum;
  ilist = comb->list->ilist;

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
