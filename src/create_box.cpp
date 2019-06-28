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

#include <cstring>
#include "create_box.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "region_prism.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal create_box command");

  if (domain->box_exist)
    error->all(FLERR,"Cannot create_box after simulation box is defined");
  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");

  domain->box_exist = 1;

  // region check

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR,"Create_box region ID does not exist");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Create_box region does not support a bounding box");

  domain->regions[iregion]->init();

  // if region not prism:
  //   setup orthogonal domain
  //   set simulation domain from region extent
  // if region is prism:
  //   seutp triclinic domain
  //   set simulation domain params from prism params

  if (strcmp(domain->regions[iregion]->style,"prism") != 0) {
    domain->triclinic = 0;
    domain->boxlo[0] = domain->regions[iregion]->extent_xlo;
    domain->boxhi[0] = domain->regions[iregion]->extent_xhi;
    domain->boxlo[1] = domain->regions[iregion]->extent_ylo;
    domain->boxhi[1] = domain->regions[iregion]->extent_yhi;
    domain->boxlo[2] = domain->regions[iregion]->extent_zlo;
    domain->boxhi[2] = domain->regions[iregion]->extent_zhi;

  } else {
    domain->triclinic = 1;
    RegPrism *region = (RegPrism *) domain->regions[iregion];
    domain->boxlo[0] = region->xlo;
    domain->boxhi[0] = region->xhi;
    domain->boxlo[1] = region->ylo;
    domain->boxhi[1] = region->yhi;
    domain->boxlo[2] = region->zlo;
    domain->boxhi[2] = region->zhi;
    domain->xy = region->xy;
    domain->xz = region->xz;
    domain->yz = region->yz;
  }

  // if molecular, zero out topology info

  if (atom->molecular) {
    atom->bond_per_atom = 0;
    atom->angle_per_atom = 0;
    atom->dihedral_per_atom = 0;
    atom->improper_per_atom = 0;
    atom->nbonds = 0;
    atom->nangles = 0;
    atom->ndihedrals = 0;
    atom->nimpropers = 0;
  }

  // set atom and topology type quantities

  atom->ntypes = force->inumeric(FLERR,arg[0]);
  atom->nbondtypes = 0;
  atom->nangletypes = 0;
  atom->ndihedraltypes = 0;
  atom->nimpropertypes = 0;

  // process optional args that can overwrite default settings

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bond/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->bonds_allow)
        error->all(FLERR,"No bonds allowed with this atom style");
      atom->nbondtypes = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"angle/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->angles_allow)
        error->all(FLERR,"No angles allowed with this atom style");
      atom->nangletypes = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dihedral/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->dihedrals_allow)
        error->all(FLERR,"No dihedrals allowed with this atom style");
      atom->ndihedraltypes = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"improper/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->impropers_allow)
        error->all(FLERR,"No impropers allowed with this atom style");
      atom->nimpropertypes = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/bond/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->bonds_allow)
        error->all(FLERR,"No bonds allowed with this atom style");
      atom->bond_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/angle/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->angles_allow)
        error->all(FLERR,"No angles allowed with this atom style");
      atom->angle_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/dihedral/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->dihedrals_allow)
        error->all(FLERR,"No dihedrals allowed with this atom style");
      atom->dihedral_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/improper/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      if (!atom->avec->impropers_allow)
        error->all(FLERR,"No impropers allowed with this atom style");
      atom->improper_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/special/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_box command");
      force->special_extra = force->inumeric(FLERR,arg[iarg+1]);
      atom->maxspecial += force->special_extra;
      iarg += 2;
    } else error->all(FLERR,"Illegal create_box command");
  }

  // problem setup using info from header
  // deallocate/grow insures any extra settings are used for topology arrays
  // necessary in case no create_atoms is performed

  update->ntimestep = 0;

  atom->allocate_type_arrays();
  atom->deallocate_topology();
  atom->avec->grow(1);

  domain->print_box("Created ");
  domain->set_initial_box();
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();
}
