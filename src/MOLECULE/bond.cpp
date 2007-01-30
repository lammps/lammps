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
#include "bond.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* -----------------------------------------------------------------------
   set bond contribution to Vdwl energy to 0.0
   a particular bond style can override this
------------------------------------------------------------------------- */

Bond::Bond(LAMMPS *lmp) : Pointers(lmp)
{
  allocated = 0;
  eng_vdwl = 0.0;
}

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void Bond::init()
{
  if (!allocated) error->all("Bond coeffs are not set");
  for (int i = 1; i <= atom->nbondtypes; i++)
    if (setflag[i] == 0) error->all("All bond coeffs are not set");
  init_style();
}
