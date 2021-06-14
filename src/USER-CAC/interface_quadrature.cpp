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

#include "interface_quadrature.h"
#include <cstring>
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

InterfaceQuadrature::InterfaceQuadrature(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void InterfaceQuadrature::command(int narg, char **arg)
{
  //check if simulation box has been defined
  if (domain->box_exist == 0)
    error->all(FLERR,"interface_quadrature command before simulation box is defined");
  //check if atom and element data has been read in already
  if (atom->natoms == 0)
    error->all(FLERR,"interface_quadrature command before material data has been defined");
  //check if CAC atom style is defined
  if(!atom->CAC_flag)
  error->all(FLERR, "interface_quadrature command requires a CAC atom style");

  if (strcmp(arg[0], "on") == 0) atom->interface_quadrature = 1;
  else if (strcmp(arg[0], "off") == 0) atom->interface_quadrature = 0;
  else error->all(FLERR, "Unexpected argument in interface_quadrature command");
}
