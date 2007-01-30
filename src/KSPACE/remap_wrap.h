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

#ifndef REMAP_WRAP_H
#define REMAP_WRAP_H

#include "pointers.h"
#include "remap.h"

namespace LAMMPS_NS {

class Remap : protected Pointers {
 public:
  Remap(class LAMMPS *, MPI_Comm,int,int,int,int,int,int,
	int,int,int,int,int,int,int,int,int,int);
  ~Remap();
  void perform(double *, double *, double *);

 private:
  struct remap_plan_3d *plan;
};

}

#endif
