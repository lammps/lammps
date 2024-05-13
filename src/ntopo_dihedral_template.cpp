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

#include "ntopo_dihedral_template.h"

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

NTopoDihedralTemplate::NTopoDihedralTemplate(LAMMPS *lmp) : NTopo(lmp)
{
  allocate_dihedral();
}

/* ---------------------------------------------------------------------- */

void NTopoDihedralTemplate::build()
{
  int i, m, atom1, atom2, atom3, atom4;
  int imol, iatom;
  tagint tagprev;
  int *num_dihedral;
  tagint **dihedral_atom1, **dihedral_atom2, **dihedral_atom3, **dihedral_atom4;
  int **dihedral_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_dihedral = onemols[imol]->num_dihedral;
    dihedral_atom1 = onemols[imol]->dihedral_atom1;
    dihedral_atom2 = onemols[imol]->dihedral_atom2;
    dihedral_atom3 = onemols[imol]->dihedral_atom3;
    dihedral_atom4 = onemols[imol]->dihedral_atom4;
    dihedral_type = onemols[imol]->dihedral_type;

    for (m = 0; m < num_dihedral[iatom]; m++) {
      atom1 = atom->map(dihedral_atom1[iatom][m] + tagprev);
      atom2 = atom->map(dihedral_atom2[iatom][m] + tagprev);
      atom3 = atom->map(dihedral_atom3[iatom][m] + tagprev);
      atom4 = atom->map(dihedral_atom4[iatom][m] + tagprev);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR)
          error->one(FLERR, "Dihedral atoms {} {} {} {} missing on proc {} at step {}",
                     dihedral_atom1[iatom][m] + tagprev, dihedral_atom2[iatom][m] + tagprev,
                     dihedral_atom3[iatom][m] + tagprev, dihedral_atom4[iatom][m] + tagprev, me,
                     update->ntimestep);
        continue;
      }
      atom1 = domain->closest_image(i, atom1);
      atom2 = domain->closest_image(i, atom2);
      atom3 = domain->closest_image(i, atom3);
      atom4 = domain->closest_image(i, atom4);
      if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
        if (ndihedrallist == maxdihedral) {
          maxdihedral += DELTA;
          memory->grow(dihedrallist, maxdihedral, 5, "neigh_topo:dihedrallist");
        }
        dihedrallist[ndihedrallist][0] = atom1;
        dihedrallist[ndihedrallist][1] = atom2;
        dihedrallist[ndihedrallist][2] = atom3;
        dihedrallist[ndihedrallist][3] = atom4;
        dihedrallist[ndihedrallist][4] = dihedral_type[iatom][m];
        ndihedrallist++;
      }
    }
  }

  if (cluster_check) dihedral_check(ndihedrallist, dihedrallist);
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing, &all, 1, MPI_INT, MPI_SUM, world);
  if (all && (me == 0))
    error->warning(FLERR, "Dihedral atoms missing at step {}", update->ntimestep);
}
