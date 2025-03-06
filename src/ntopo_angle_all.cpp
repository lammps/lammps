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

#include "ntopo_angle_all.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "output.h"
#include "thermo.h"
#include "update.h"

using namespace LAMMPS_NS;

static constexpr int DELTA = 10000;

/* ---------------------------------------------------------------------- */

NTopoAngleAll::NTopoAngleAll(LAMMPS *lmp) : NTopo(lmp)
{
  allocate_angle();
}

/* ---------------------------------------------------------------------- */

void NTopoAngleAll::build()
{
  int i, m, atom1, atom2, atom3;

  int nlocal = atom->nlocal;
  int *num_angle = atom->num_angle;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nanglelist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_angle[i]; m++) {
      atom1 = atom->map(angle_atom1[i][m]);
      atom2 = atom->map(angle_atom2[i][m]);
      atom3 = atom->map(angle_atom3[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR)
          error->one(FLERR, "Angle atoms {} {} {} missing on proc {} at step {}", angle_atom1[i][m],
                     angle_atom2[i][m], angle_atom3[i][m], me, update->ntimestep);
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
        anglelist[nanglelist][3] = angle_type[i][m];
        nanglelist++;
      }
    }

  if (cluster_check) angle_check();
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing, &all, 1, MPI_INT, MPI_SUM, world);
  if (all && (me == 0)) error->warning(FLERR, "Angle atoms missing at step {}", update->ntimestep);
}
