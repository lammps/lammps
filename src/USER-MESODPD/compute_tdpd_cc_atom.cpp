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

#include <cstring>
#include "compute_tdpd_cc_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTDPDCCAtom::ComputeTDPDCCAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Number of arguments for compute tdpd/cc/atom command != 4");
  if (atom->tdpd_flag != 1) error->all(FLERR,"compute tdpd/cc/atom command requires atom_style with concentration (e.g. tdpd)");

  index = force->inumeric(FLERR,arg[3]);

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  cc_vector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeTDPDCCAtom::~ComputeTDPDCCAtom()
{
  memory->sfree(cc_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeTDPDCCAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"cc_vector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute cc_vector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeTDPDCCAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow cc_vector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(cc_vector);
    nmax = atom->nmax;
    cc_vector = (double *) memory->smalloc(nmax*sizeof(double),"cc_vector/atom:cc_vector");
    vector_atom = cc_vector;
  }

  double **cc = atom->cc;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
         cc_vector[i] = cc[i][index-1];
      }
      else
         cc_vector[i] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeTDPDCCAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
