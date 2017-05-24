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
   			 Oleg Sergeev (VNIIA, sergeev@vniia.ru)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "stdlib.h"
#include "math.h"
#include "atom.h"
#include "string.h"
#include "fix_ave_atom.h"
#include "fix_reaxc_species_omp.h"
#include "domain.h"
#include "update.h"
#include "pair_reaxc_omp.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "reaxc_list.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCSpeciesOMP::FixReaxCSpeciesOMP(LAMMPS *lmp, int narg, char **arg) :
  FixReaxCSpecies(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpeciesOMP::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix reax/c/species unless atoms have IDs");

  reaxc = (PairReaxCOMP *) force->pair_match("reax/c/omp",1);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/species/omp without "
		  "pair_style reax/c/omp");

  reaxc->fixspecies_flag = 1;
  nvalid = update->ntimestep+nfreq;

  // check if this fix has been called twice
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"reax/c/species/omp") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one fix reax/c/species");

  if (!setupflag) {
    // create a compute to store properties
    create_compute();

    // create a fix to point to fix_ave_atom for averaging stored properties
    create_fix();

    setupflag = 1;
  }

}
