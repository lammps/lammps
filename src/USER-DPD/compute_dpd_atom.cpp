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
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include "compute_dpd_atom.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDpdAtom::ComputeDpdAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), dpdAtom(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute dpd/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;

  nmax = 0;

  if (atom->dpd_flag != 1) error->all(FLERR,"compute dpd requires atom_style with internal temperature and energies (e.g. dpd)");
}

/* ---------------------------------------------------------------------- */

ComputeDpdAtom::~ComputeDpdAtom()
{
  memory->destroy(dpdAtom);
}

/* ---------------------------------------------------------------------- */

void ComputeDpdAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"dpd/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute dpd/atom command");
}

/* ----------------------------------------------------------------------
   gather compute vector data from other nodes
------------------------------------------------------------------------- */

void ComputeDpdAtom::compute_peratom()
{

  invoked_peratom = update->ntimestep;

  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *uChem = atom->uChem;
  double *dpdTheta = atom->dpdTheta;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  if (nlocal > nmax) {
    memory->destroy(dpdAtom);
    nmax = atom->nmax;
    memory->create(dpdAtom,nmax,size_peratom_cols,"dpd/atom:dpdAtom");
    array_atom = dpdAtom;
  }

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit){
      dpdAtom[i][0] =  uCond[i];
      dpdAtom[i][1] =  uMech[i];
      dpdAtom[i][2] =  uChem[i];
      dpdAtom[i][3] =  dpdTheta[i];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDpdAtom::memory_usage()
{
  double bytes = size_peratom_cols * nmax * sizeof(double);
  return bytes;
}
