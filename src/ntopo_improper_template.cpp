/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Templatel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "ntopo_improper_template.h"
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

NTopoImproperTemplate::NTopoImproperTemplate(LAMMPS *lmp) :
  NTopo(lmp)
{
  allocate_improper();
}

/* ---------------------------------------------------------------------- */

void NTopoImproperTemplate::build()
{
  int i,m,atom1,atom2,atom3,atom4;
  int imol,iatom;
  tagint tagprev;
  int *num_improper;
  tagint **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;
  int **improper_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nimproperlist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_improper = onemols[imol]->num_improper;
    improper_atom1 = onemols[imol]->improper_atom1;
    improper_atom2 = onemols[imol]->improper_atom2;
    improper_atom3 = onemols[imol]->improper_atom3;
    improper_atom4 = onemols[imol]->improper_atom4;
    improper_type = onemols[imol]->improper_type;

    for (m = 0; m < num_improper[iatom]; m++) {
      atom1 = atom->map(improper_atom1[iatom][m]+tagprev);
      atom2 = atom->map(improper_atom2[iatom][m]+tagprev);
      atom3 = atom->map(improper_atom3[iatom][m]+tagprev);
      atom4 = atom->map(improper_atom4[iatom][m]+tagprev);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == Thermo::ERROR) {
          char str[128];
          sprintf(str,"Improper atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " "
                  TAGINT_FORMAT " " TAGINT_FORMAT
                  " missing on proc %d at step " BIGINT_FORMAT,
                  improper_atom1[iatom][m]+tagprev,
                  improper_atom2[iatom][m]+tagprev,
                  improper_atom3[iatom][m]+tagprev,
                  improper_atom4[iatom][m]+tagprev,
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain-> closest_image(i,atom4);
      if (newton_bond ||
          (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
        if (nimproperlist == maximproper) {
          maximproper += DELTA;
          memory->grow(improperlist,maximproper,5,"neigh_topo:improperlist");
        }
        improperlist[nimproperlist][0] = atom1;
        improperlist[nimproperlist][1] = atom2;
        improperlist[nimproperlist][2] = atom3;
        improperlist[nimproperlist][3] = atom4;
        improperlist[nimproperlist][4] = improper_type[iatom][m];
        nimproperlist++;
      }
    }
  }

  if (cluster_check) dihedral_check(nimproperlist,improperlist);
  if (lostbond == Thermo::IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Improper atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}
