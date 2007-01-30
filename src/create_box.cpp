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

#include "stdlib.h"
#include "string.h"
#include "create_box.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal create_box command");

  if (domain->box_exist) 
    error->all("Cannot create_box after simulation box is defined");
  if (force->dimension == 2 && domain->zperiodic == 0)
    error->all("Cannot run 2d simulation with nonperiodic Z dimension");

  // find region ID

  int iregion;
  for (iregion = 0; iregion < domain->nregion; iregion++)
    if (strcmp(arg[1],domain->regions[iregion]->id) == 0) break;
  if (iregion == domain->nregion) 
    error->all("Create_box region ID does not exist");

  if (domain->regions[iregion]->interior == 0)
    error->all("Create_box region must be of type inside");

  // set global box from region extent

  domain->boxxlo = domain->regions[iregion]->extent_xlo;
  domain->boxxhi = domain->regions[iregion]->extent_xhi;
  domain->boxylo = domain->regions[iregion]->extent_ylo;
  domain->boxyhi = domain->regions[iregion]->extent_yhi;
  domain->boxzlo = domain->regions[iregion]->extent_zlo;
  domain->boxzhi = domain->regions[iregion]->extent_zhi;
  
  domain->box_exist = 1;

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Created box = (%g %g %g) to (%g %g %g)\n",
	      domain->boxxlo,domain->boxylo,domain->boxzlo,
	      domain->boxxhi,domain->boxyhi,domain->boxzhi);
    if (logfile)
      fprintf(logfile,"Created box = (%g %g %g) to (%g %g %g)\n",
	      domain->boxxlo,domain->boxylo,domain->boxzlo,
	      domain->boxxhi,domain->boxyhi,domain->boxzhi);
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

  atom->ntypes = atoi(arg[0]);
  atom->nbondtypes = 0;
  atom->nangletypes = 0;
  atom->ndihedraltypes = 0;
  atom->nimpropertypes = 0;

  // problem setup using info from header
  // no call to atom->grow since create_atoms or fixes will do it

  update->ntimestep = 0;

  atom->allocate_type_arrays();

  domain->set_initial_box();
  domain->set_global_box();
  comm->set_procs();
  domain->set_local_box();
}
