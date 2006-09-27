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
#include "angle.h"
#include "atom.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Angle::Angle()
{
  allocated = 0;
  PI = 4.0*atan(1.0);
}

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void Angle::init()
{
  if (!allocated) error->all("Angle coeffs are not set");
  for (int i = 1; i <= atom->nangletypes; i++)
    if (setflag[i] == 0) error->all("All angle coeffs are not set");
}

