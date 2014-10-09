/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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

#ifdef COMPUTE_CLASS

ComputeStyle(temp/cuda,ComputeTempCuda)

#else

#ifndef LMP_COMPUTE_TEMP_CUDA_H
#define LMP_COMPUTE_TEMP_CUDA_H

#include "compute.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class ComputeTempCuda : public Compute {
 public:
  ComputeTempCuda(class LAMMPS *, int, char **);
  ~ComputeTempCuda();
  void init() {}
  void setup();
  double compute_scalar();
  void compute_vector();

 private:
  class Cuda *cuda;
  int fix_dof;
  double tfactor;

  void dof_compute();
  double t_vector[6];
  double t_scalar;
  cCudaData<double     , ENERGY_CFLOAT                   , x>* cu_t_scalar;
  cCudaData<double     , ENERGY_CFLOAT                   , x>* cu_t_vector;

};

}

#endif
#endif
