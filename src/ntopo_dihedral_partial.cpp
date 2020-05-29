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

#include "ntopo_dihedral_partial.h"
#include <mpi.h>
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

NTopoDihedralPartial::NTopoDihedralPartial(LAMMPS *lmp) :
  NTopo(lmp)
{
  allocate_dihedral();
}

/* ---------------------------------------------------------------------- */

void NTopoDihedralPartial::build()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_dihedral[i]; m++) {
      if (dihedral_type[i][m] <= 0) continue;
      atom1 = atom->map(dihedral_atom1[i][m]);
      atom2 = atom->map(dihedral_atom2[i][m]);
      atom3 = atom->map(dihedral_atom3[i][m]);
      atom4 = atom->map(dihedral_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR) {
          char str[128];
          sprintf(str,"Dihedral atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " "
                  TAGINT_FORMAT " " TAGINT_FORMAT
                  " missing on proc %d at step " BIGINT_FORMAT,
                  dihedral_atom1[i][m],dihedral_atom2[i][m],
                  dihedral_atom3[i][m],dihedral_atom4[i][m],
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
      if (newton_bond ||
          (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
        if (ndihedrallist == maxdihedral) {
          maxdihedral += DELTA;
          memory->grow(dihedrallist,maxdihedral,5,"neigh_topo:dihedrallist");
        }
        dihedrallist[ndihedrallist][0] = atom1;
        dihedrallist[ndihedrallist][1] = atom2;
        dihedrallist[ndihedrallist][2] = atom3;
        dihedrallist[ndihedrallist][3] = atom4;
        dihedrallist[ndihedrallist][4] = dihedral_type[i][m];
        ndihedrallist++;
      }
    }

  if (cluster_check) dihedral_check(ndihedrallist,dihedrallist);
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Dihedral atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}
