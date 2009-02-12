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

#ifndef COMPUTE_COORD_ATOM_H
#define COMPUTE_COORD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCoordAtom : public Compute {
 public:
  ComputeCoordAtom(class LAMMPS *, int, char **);
  ~ComputeCoordAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double cutsq;
  class NeighList *list;
  double *coordination;
};

}

#endif
