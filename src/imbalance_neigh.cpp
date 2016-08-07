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


#include "pointers.h"
#include "imbalance_neigh.h"
#include "atom.h"
#include "error.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

int ImbalanceNeigh::options(int narg, char **arg)
{
  Error *error = _lmp->error;
  Force *force = _lmp->force;

  if (narg < 1) error->all(FLERR,"Illegal balance weight command");
  _factor = force->numeric(FLERR,arg[0]);
  if (_factor < 0.0) error->all(FLERR,"Illegal balance weight command");
  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceNeigh::compute(double *weight)
{
  // find suitable existing neighbor list
  Error *error = _lmp->error;
  Neighbor *neighbor = _lmp->neighbor;
  int req;

  // we can only make use of certain (conventional) neighbor lists.
  for (req = 0; req < neighbor->old_nrequest; ++req) {
    if ((neighbor->old_requests[req]->half ||
         neighbor->old_requests[req]->gran ||
         neighbor->old_requests[req]->respaouter ||
         neighbor->old_requests[req]->half_from_full) &&
        neighbor->old_requests[req]->skip == 0 &&
        neighbor->lists[req] && neighbor->lists[req]->numneigh) break;
  }

  if (req >= neighbor->old_nrequest || neighbor->ago < 0) {
    if (_lmp->comm->me == 0 && !did_warn)
      error->warning(FLERR,"No suitable neighbor list found. "
                     "Neighbor weighted balancing skipped");
    did_warn = 1;
    return;
  }

  if (_factor > 0.0) {

    MPI_Comm world = _lmp->world;
    NeighList *list = neighbor->lists[req];
    bigint neighsum = 0;

    const int inum = list->inum;
    const int * const ilist = list->ilist;
    const int * const numneigh = list->numneigh;

    // first pass: get local number of neighbors.
    for (int i = 0; i < inum; ++i) neighsum += numneigh[ilist[i]];

    double allatoms = static_cast <double>(_lmp->atom->natoms);
    if (allatoms == 0.0) allatoms = 1.0;
    double allavg;
    double myavg = static_cast<double>(neighsum)/allatoms;
    MPI_Allreduce(&myavg,&allavg,1,MPI_DOUBLE,MPI_SUM,world);

    // second pass: compute and apply weights
    double scale = 1.0/allavg;
    for (int ii = 0; ii < inum; ++ii) {
      const int i = ilist[ii];
      weight[i] *= (1.0-_factor) + _factor*scale*numneigh[i];
    }
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceNeigh::info(FILE *fp)
{
  if (_factor > 0.0)
    fprintf(fp,"  neigh weight factor: %g\n",_factor);
}
