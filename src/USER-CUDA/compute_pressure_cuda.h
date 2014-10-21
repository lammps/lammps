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

ComputeStyle(pressure/cuda,ComputePressureCuda)

#else

#ifndef LMP_COMPUTE_PRESSURE_CUDA_H
#define LMP_COMPUTE_PRESSURE_CUDA_H

#include "compute_pressure.h"

namespace LAMMPS_NS {

class ComputePressureCuda : public ComputePressure {
 public:
  ComputePressureCuda(class LAMMPS *, int, char **);
  ~ComputePressureCuda() {}
  double compute_scalar();
  void compute_vector();

  private:
  class Cuda *cuda;
};

}

#endif
#endif
