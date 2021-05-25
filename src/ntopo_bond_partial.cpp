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

#include "ntopo_bond_partial.h"

#include "atom.h"
#include "force.h"
#include "domain.h"
#include "update.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

NTopoBondPartial::NTopoBondPartial(LAMMPS *lmp) : NTopo(lmp)
{
  allocate_bond();
}

/* ---------------------------------------------------------------------- */

void NTopoBondPartial::build()
{
  int i,m,atom1;

  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nbondlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_type[i][m] <= 0) continue;
      atom1 = atom->map(bond_atom[i][m]);
      if (atom1 == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR)
          error->one(FLERR,"Bond atoms {} {} missing on "
                                       "proc {} at step {}",tag[i],
                                       bond_atom[i][m],me,update->ntimestep);
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      if (newton_bond || i < atom1) {
        if (nbondlist == maxbond) {
          maxbond += DELTA;
          memory->grow(bondlist,maxbond,3,"neigh_topo:bondlist");
        }
        bondlist[nbondlist][0] = i;
        bondlist[nbondlist][1] = atom1;
        bondlist[nbondlist][2] = bond_type[i][m];
        nbondlist++;
      }
    }

  if (cluster_check) bond_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all && (me == 0))
    error->warning(FLERR,"Bond atoms missing at step {}",update->ntimestep);
}
