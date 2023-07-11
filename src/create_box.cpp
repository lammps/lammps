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

#include "create_box.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "region.h"
#include "region_prism.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "create_box", error);

  if (domain->box_exist) error->all(FLERR, "Cannot create_box after simulation box is defined");
  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR, "Cannot run 2d simulation with nonperiodic Z dimension");

  domain->box_exist = 1;

  // region check

  auto region = domain->get_region_by_id(arg[1]);
  if (!region) error->all(FLERR, "Create_box region {} does not exist", arg[1]);
  if (region->bboxflag == 0) error->all(FLERR, "Create_box region does not support a bounding box");

  region->init();

  // if region not prism:
  //   setup orthogonal domain
  //   set simulation domain from region extent
  // if region is prism:
  //   seutp triclinic domain
  //   set simulation domain params from prism params

  if (strcmp(region->style, "prism") != 0) {
    domain->triclinic = 0;
    domain->boxlo[0] = region->extent_xlo;
    domain->boxhi[0] = region->extent_xhi;
    domain->boxlo[1] = region->extent_ylo;
    domain->boxhi[1] = region->extent_yhi;
    domain->boxlo[2] = region->extent_zlo;
    domain->boxhi[2] = region->extent_zhi;

  } else {
    domain->triclinic = 1;
    auto prism = dynamic_cast<RegPrism *>(region);
    domain->boxlo[0] = prism->xlo;
    domain->boxhi[0] = prism->xhi;
    domain->boxlo[1] = prism->ylo;
    domain->boxhi[1] = prism->yhi;
    domain->boxlo[2] = prism->zlo;
    domain->boxhi[2] = prism->zhi;
    domain->xy = prism->xy;
    domain->xz = prism->xz;
    domain->yz = prism->yz;
  }

  // if molecular, zero out topology info

  if (atom->molecular != Atom::ATOMIC) {
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

  atom->ntypes = utils::inumeric(FLERR, arg[0], false, lmp);
  atom->nbondtypes = 0;
  atom->nangletypes = 0;
  atom->ndihedraltypes = 0;
  atom->nimpropertypes = 0;

  // process optional args that can overwrite default settings

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "bond/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_box bond/type", error);
      if (!atom->avec->bonds_allow)
        error->all(FLERR, "No bonds allowed with atom style {}", atom->get_style());
      atom->nbondtypes = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "angle/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_box angle/types", error);
      if (!atom->avec->angles_allow)
        error->all(FLERR, "No angles allowed with atom style {}", atom->get_style());
      atom->nangletypes = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dihedral/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_box dihedral/types", error);
      if (!atom->avec->dihedrals_allow)
        error->all(FLERR, "No dihedrals allowed with atom style {}", atom->get_style());
      atom->ndihedraltypes = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "improper/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_box improper/types", error);
      if (!atom->avec->impropers_allow)
        error->all(FLERR, "No impropers allowed with atom style {}", atom->get_style());
      atom->nimpropertypes = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/bond/per/atom") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_box extra/bond/per/atom", error);
      if (!atom->avec->bonds_allow)
        error->all(FLERR, "No bonds allowed with atom style {}", atom->get_style());
      atom->bond_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/angle/per/atom") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_box extra/angle/per/atom", error);
      if (!atom->avec->angles_allow)
        error->all(FLERR, "No angles allowed with atom style {}", atom->get_style());
      atom->angle_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/dihedral/per/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "create_box extra/dihedral/per/atom", error);
      if (!atom->avec->dihedrals_allow)
        error->all(FLERR, "No dihedrals allowed with atom style {}", atom->get_style());
      atom->dihedral_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/improper/per/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "create_box extra/improper/per/atom", error);
      if (!atom->avec->impropers_allow)
        error->all(FLERR, "No impropers allowed with atom style {}", atom->get_style());
      atom->improper_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/special/per/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "create_box extra/special/per/atom", error);
      force->special_extra = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      atom->maxspecial += force->special_extra;
      iarg += 2;
    } else
      error->all(FLERR, "Unknown create_box keyword: {}", arg[iarg]);
  }

  // problem setup using info from header
  // deallocate/grow ensures any extra settings are used for topology arrays
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
