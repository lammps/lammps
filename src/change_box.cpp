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

#include "string.h"
#include "change_box.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "output.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{ORTHO,TRICLINIC};

/* ---------------------------------------------------------------------- */

ChangeBox::ChangeBox(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void ChangeBox::command(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("Change_box command before simulation box is defined");
  if (narg != 1) error->all("Illegal change_box command");

  int style;
  if (strcmp(arg[0],"ortho") == 0) style = ORTHO;
  else if (strcmp(arg[0],"triclinic") == 0) style = TRICLINIC;

  if (style == ORTHO && domain->triclinic == 0)
    error->all("Change_box operation is invalid");
  if (style == TRICLINIC && domain->triclinic == 1)
    error->all("Change_box operation is invalid");
  if (style == ORTHO && 
      (domain->xy != 0.0 || domain->yz != 0.0 || domain->xz != 0.0))
      error->all("Cannot change box to orthogonal when tilt is non-zero");

  if (output->ndump)
    error->all("Cannot change box with dumps defined");
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->no_change_box)
      error->all("Cannot change box with certain fixes defined");

  if (style == ORTHO) domain->triclinic = 0;
  else domain->triclinic = 1;
  domain->xy = domain->yz = domain->xz = 0.0;

  domain->set_global_box();
  if (style == TRICLINIC) domain->set_lamda_box();
  domain->set_local_box();
}
