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

#include "math.h"
#include "dihedral.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   set dihedral contribution to Vdwl and Coulombic energy to 0.0
   DihedralCharmm will override this
------------------------------------------------------------------------- */

Dihedral::Dihedral(LAMMPS *lmp) : Pointers(lmp)
{
  allocated = 0;
  eng_vdwl = eng_coul = 0.0;
  PI = 4.0*atan(1.0);
}

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void Dihedral::init()
{
  if (!allocated) error->all("Dihedral coeffs are not set");
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    if (setflag[i] == 0) error->all("All dihedral coeffs are not set");
  init_style();
}

