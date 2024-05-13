/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "ntopo_angle_template.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "molecule.h"
#include "output.h"
#include "thermo.h"
#include "update.h"

using namespace LAMMPS_NS;

static constexpr int DELTA = 10000;

/* ---------------------------------------------------------------------- */

NTopoAngleTemplate::NTopoAngleTemplate(LAMMPS *lmp) : NTopo(lmp)
{
  allocate_angle();
}

/* ---------------------------------------------------------------------- */

void NTopoAngleTemplate::build()
{
  int i, m, atom1, atom2, atom3;
  int imol, iatom;
  tagint tagprev;
  int *num_angle;
  tagint **angle_atom1, **angle_atom2, **angle_atom3;
  int **angle_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nanglelist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_angle = onemols[imol]->num_angle;
    angle_atom1 = onemols[imol]->angle_atom1;
    angle_atom2 = onemols[imol]->angle_atom2;
    angle_atom3 = onemols[imol]->angle_atom3;
    angle_type = onemols[imol]->angle_type;

    for (m = 0; m < num_angle[iatom]; m++) {
      if (angle_type[iatom][m] <= 0) continue;
      atom1 = atom->map(angle_atom1[iatom][m] + tagprev);
      atom2 = atom->map(angle_atom2[iatom][m] + tagprev);
      atom3 = atom->map(angle_atom3[iatom][m] + tagprev);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR)
          error->one(FLERR, "Angle atoms {} {} {} missing on proc {} at step {}",
                     angle_atom1[iatom][m] + tagprev, angle_atom2[iatom][m] + tagprev,
                     angle_atom3[iatom][m] + tagprev, me, update->ntimestep);
        continue;
      }
      atom1 = domain->closest_image(i, atom1);
      atom2 = domain->closest_image(i, atom2);
      atom3 = domain->closest_image(i, atom3);
      if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
        if (nanglelist == maxangle) {
          maxangle += DELTA;
          memory->grow(anglelist, maxangle, 4, "neigh_topo:anglelist");
        }
        anglelist[nanglelist][0] = atom1;
        anglelist[nanglelist][1] = atom2;
        anglelist[nanglelist][2] = atom3;
        anglelist[nanglelist][3] = angle_type[iatom][m];
        nanglelist++;
      }
    }
  }

  if (cluster_check) angle_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing, &all, 1, MPI_INT, MPI_SUM, world);
  if (all && (me == 0)) error->warning(FLERR, "Angle atoms missing at step {}", update->ntimestep);
}
