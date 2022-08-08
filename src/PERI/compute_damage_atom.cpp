// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Parks (SNL)
------------------------------------------------------------------------- */

#include "compute_damage_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_peri_neigh.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDamageAtom::ComputeDamageAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), damage(nullptr), fix_peri_neigh(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal compute damage/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeDamageAtom::~ComputeDamageAtom()
{
  memory->destroy(damage);
}

/* ---------------------------------------------------------------------- */

void ComputeDamageAtom::init()
{
  if ((comm->me == 0) && (modify->get_compute_by_style("damage/atom").size() > 1))
    error->warning(FLERR,"More than one compute dilatation/atom");

  // find associated PERI_NEIGH fix that must exist

  auto fixes = modify->get_fix_by_style("PERI_NEIGH");
  if (fixes.size() == 0)
    error->all(FLERR,"Compute damage/atom requires a peridynamic potential");
  else fix_peri_neigh = dynamic_cast<FixPeriNeigh *>(fixes.front());
}

/* ---------------------------------------------------------------------- */

void ComputeDamageAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow damage array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(damage);
    nmax = atom->nmax;
    memory->create(damage,nmax,"damage/atom:damage");
    vector_atom = damage;
  }

  // compute damage for each atom in group

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *vfrac = atom->vfrac;
  double *vinter = fix_peri_neigh->vinter;
  tagint **partner = fix_peri_neigh->partner;
  int *npartner = fix_peri_neigh->npartner;
  int i,j,jj,jnum;

  double damage_temp;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      jnum = npartner[i];
      damage_temp = 0.0;
      for (jj = 0; jj < jnum; jj++) {
        if (partner[i][jj] == 0) continue;

        // look up local index of this partner particle
        // skip if particle is "lost"

        j = atom->map(partner[i][jj]);
        if (j < 0) continue;

        damage_temp += vfrac[j];
      }

      if (vinter[i] != 0.0) damage[i] = 1.0 - damage_temp/vinter[i];
      else damage[i] = 0.0;

    } else damage[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDamageAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
