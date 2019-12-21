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
   Contributing author: Rezwanur Rahman, John Foster (UTSA)
------------------------------------------------------------------------- */

#include "compute_dilatation_atom.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "fix.h"
#include "force.h"
#include "pair_peri_lps.h"
#include "pair_peri_ves.h"
#include "pair_peri_eps.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDilatationAtom::
ComputeDilatationAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute Dilatation/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  dilatation = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDilatationAtom::~ComputeDilatationAtom()
{
  memory->destroy(dilatation);
}

/* ---------------------------------------------------------------------- */

void ComputeDilatationAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"dilatation/peri") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute dilatation/atom");

  // check PD pair style

  isPMB = isLPS = isVES = isEPS = 0;
  if (force->pair_match("^peri/pmb",0)) isPMB = 1;
  if (force->pair_match("^peri/lps",0)) isLPS = 1;
  if (force->pair_match("^peri/ves",0)) isVES = 1;
  if (force->pair_match("^peri/eps",0)) isEPS = 1;

  if (isPMB)
    error->all(FLERR,"Compute dilatation/atom cannot be used "
               "with this pair style");

  // find associated PERI_NEIGH fix that must exist

  if (modify->find_fix_by_style("^PERI_NEIGH") == -1)
    error->all(FLERR,"Compute dilatation/atom requires Peridynamic pair style");
}

/* ---------------------------------------------------------------------- */

void ComputeDilatationAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow dilatation array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(dilatation);
    nmax = atom->nmax;
    memory->create(dilatation,nmax,"dilatation/atom:dilatation");
    vector_atom = dilatation;
  }

  // extract dilatation for each atom in group

  double *theta;
  Pair *anypair = force->pair_match("peri",0);
  if (isLPS) theta = ((PairPeriLPS *) anypair)->theta;
  if (isVES) theta = ((PairPeriVES *) anypair)->theta;
  if (isEPS) theta = ((PairPeriEPS *) anypair)->theta;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dilatation[i] = theta[i];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDilatationAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
