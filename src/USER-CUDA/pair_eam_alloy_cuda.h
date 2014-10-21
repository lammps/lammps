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

#ifdef PAIR_CLASS

PairStyle(eam/alloy/cuda,PairEAMAlloyCuda)

#else

#ifndef LMP_PAIR_EAM_CUDA_ALLOY_H
#define LMP_PAIR_EAM_CUDA_ALLOY_H

#include "pair_eam_cuda.h"

namespace LAMMPS_NS {

// use virtual public since this class is parent in multiple inheritance

class PairEAMAlloyCuda : virtual public PairEAMCuda {
 public:
  PairEAMAlloyCuda(class LAMMPS *);
  virtual ~PairEAMAlloyCuda() {}
  void coeff(int, char **);

 protected:
  class Cuda *cuda;
  void read_file(char *);
  void file2array();
};

}

#endif
#endif
