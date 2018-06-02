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

#include <mpi.h>
#include "imbalance_neigh.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* -------------------------------------------------------------------- */

ImbalanceNeigh::ImbalanceNeigh(LAMMPS *lmp) : Imbalance(lmp)
{
  did_warn = 0;
}

/* -------------------------------------------------------------------- */

int ImbalanceNeigh::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal balance weight command");
  factor = force->numeric(FLERR,arg[0]);
  if (factor <= 0.0) error->all(FLERR,"Illegal balance weight command");
  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceNeigh::compute(double *weight)
{
  int req;

  if (factor == 0.0) return;

  // find suitable neighbor list
  // can only use certain conventional neighbor lists
  // NOTE: why not full list, if half does not exist?

  for (req = 0; req < neighbor->old_nrequest; ++req) {
    if (neighbor->old_requests[req]->half &&
        neighbor->old_requests[req]->skip == 0 &&
        neighbor->lists[req] && neighbor->lists[req]->numneigh) break;
  }

  if (req >= neighbor->old_nrequest || neighbor->ago < 0) {
    if (comm->me == 0 && !did_warn)
      error->warning(FLERR,"Balance weight neigh skipped b/c no list found");
    did_warn = 1;
    return;
  }

  // neighsum = total neigh count for atoms on this proc
  // localwt = weight assigned to each owned atom

  NeighList *list = neighbor->lists[req];
  const int inum = list->inum;
  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  int nlocal = atom->nlocal;

  bigint neighsum = 0;
  for (int i = 0; i < inum; ++i) neighsum += numneigh[ilist[i]];
  double localwt = 0.0;
  if (nlocal) localwt = 1.0*neighsum/nlocal;

  if (nlocal && localwt <= 0.0) error->one(FLERR,"Balance weight <= 0.0");

  // apply factor if specified != 1.0
  // wtlo,wthi = lo/hi values excluding 0.0 due to no atoms on this proc
  // lo value does not change
  // newhi = new hi value to give hi/lo ratio factor times larger/smaller
  // expand/contract all localwt values from lo->hi to lo->newhi

  if (factor != 1.0) {
    double wtlo,wthi;
    if (localwt == 0.0) localwt = BIG;
    MPI_Allreduce(&localwt,&wtlo,1,MPI_DOUBLE,MPI_MIN,world);
    if (localwt == BIG) localwt = 0.0;
    MPI_Allreduce(&localwt,&wthi,1,MPI_DOUBLE,MPI_MAX,world);
    if (wtlo == wthi) return;

    double newhi = wthi*factor;
    localwt = wtlo + ((localwt-wtlo)/(wthi-wtlo)) * (newhi-wtlo);
  }

  for (int i = 0; i < nlocal; i++) weight[i] *= localwt;
}

/* -------------------------------------------------------------------- */

void ImbalanceNeigh::info(FILE *fp)
{
  fprintf(fp,"  neigh weight factor: %g\n",factor);
}
