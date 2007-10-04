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

#ifndef COMPUTE_EPAIR_ATOM_H
#define COMPUTE_EPAIR_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEpairAtom : public Compute {
 public:
  ComputeEpairAtom(class LAMMPS *, int, char **);
  ~ComputeEpairAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int nmax,eamstyle;
  class NeighList *list;
  double *energy;
};

}

#endif
