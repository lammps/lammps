/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "improper.h"
#include "atom.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Improper::Improper()
{
  allocated = 0;
  PI = 4.0*atan(1.0);
}

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void Improper::init()
{
  if (!allocated) error->all("Improper coeffs are not set");
  for (int i = 1; i <= atom->nimpropertypes; i++)
    if (setflag[i] == 0) error->all("All improper coeffs are not set");
}
