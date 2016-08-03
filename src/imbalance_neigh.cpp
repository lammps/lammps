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

using namespace LAMMPS_NS;

int ImbalanceNeigh::options(int narg, char **arg)
{
  Error *error = _lmp->error;
  Force *force = _lmp->force;

  if (narg < 1) error->all(FLERR,"Illegal balance weight command");
  _factor = force->numeric(FLERR,arg[0]);
  if (_factor < 0.0 || _factor > 1.0)
    error->all(FLERR,"Illegal balance weight command");
  return 1;
}

/* -------------------------------------------------------------------- */
 
void ImbalanceNeigh::compute(double *weight)
{
  const int nlocal = _lmp->atom->nlocal;
  MPI_Comm world = _lmp->world;

  if (_factor > 0.0) {
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceNeigh::info(FILE *fp)
{
  if (_factor > 0.0)
    fprintf(fp,"  neigh weight factor: %g\n",_factor);
}
