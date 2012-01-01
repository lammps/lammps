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

#ifndef LMP_NEIGHBOR_CUDA_H
#define LMP_NEIGHBOR_CUDA_H

#include "neighbor.h"

namespace LAMMPS_NS {

class NeighborCuda : public Neighbor {
 public:
  NeighborCuda(class LAMMPS *);
  void init();
  int check_distance();
  void build();

 private:
  class Cuda *cuda;

  void choose_build(int, class NeighRequest *);
  typedef void (NeighborCuda::*PairPtr)(class NeighList *);
  void full_nsq_cuda(class NeighList *);
  void full_bin_cuda(class NeighList *);
};

}

#endif
