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

PairStyle(eam/fs/cuda,PairEAMFSCuda)

#else

#ifndef LMP_PAIR_EAM_FS_CUDA_H
#define LMP_PAIR_EAM_FS_CUDA_H

#include "pair_eam_cuda.h"

namespace LAMMPS_NS {

// use virtual public since this class is parent in multiple inheritance

class PairEAMFSCuda : virtual public PairEAMCuda {
 public:
  PairEAMFSCuda(class LAMMPS *);
  virtual ~PairEAMFSCuda() {}
  void coeff(int, char **);

 protected:
  class Cuda *cuda;
  void read_file(char *);
  void file2array();
};

}

#endif
#endif
