/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/omp,FixWallGranOMP)

#else

#ifndef LMP_FIX_WALL_GRAN_OMP_H
#define LMP_FIX_WALL_GRAN_OMP_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallGranOMP : public FixWallGran {

 public:
  FixWallGranOMP(class LAMMPS *, int, char **);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
};

}

#endif
#endif
