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

#ifndef REMAP_WRAP_H
#define REMAP_WRAP_H

#include "lammps.h"
#include "remap.h"

class Remap : public LAMMPS {
 public:
  Remap(MPI_Comm,int,int,int,int,int,int,
	int,int,int,int,int,int,int,int,int,int);
  ~Remap();
  void perform(double *, double *, double *);

 private:
  struct remap_plan_3d *plan;
};

#endif
