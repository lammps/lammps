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

#include "ntopo_bond_template.h"
#include <mpi.h>
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "domain.h"
#include "update.h"
#include "output.h"
#include "thermo.h"
#include "molecule.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

NTopoBondTemplate::NTopoBondTemplate(LAMMPS *lmp) : NTopo(lmp)
{
  allocate_bond();
}

/* ---------------------------------------------------------------------- */

void NTopoBondTemplate::build()
{
  int i,m,atom1;
  int imol,iatom;
  tagint tagprev;
  int *num_bond;
  tagint **bond_atom;
  int **bond_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nbondlist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_bond = onemols[imol]->num_bond;
    bond_atom = onemols[imol]->bond_atom;
    bond_type = onemols[imol]->bond_type;

    for (m = 0; m < num_bond[iatom]; m++) {
      if (bond_type[iatom][m] <= 0) continue;
      atom1 = atom->map(bond_atom[iatom][m]+tagprev);
      if (atom1 == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR) {
          char str[128];
          sprintf(str,"Bond atoms " TAGINT_FORMAT " " TAGINT_FORMAT
                  " missing on proc %d at step " BIGINT_FORMAT,
                  tag[i],bond_atom[iatom][m]+tagprev,me,update->ntimestep);
          error->one(FLERR,str);
        }
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
        bondlist[nbondlist][2] = bond_type[iatom][m];
        nbondlist++;
      }
    }
  }

  if (cluster_check) bond_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Bond atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}
