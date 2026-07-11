/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(gaussian/grid/local,ComputeGaussianGridLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_GAUSSIAN_GRID_LOCAL_H
#define LMP_COMPUTE_GAUSSIAN_GRID_LOCAL_H

#include "compute_grid_local.h"

namespace LAMMPS_NS {

class ComputeGaussianGridLocal : public ComputeGridLocal {
 public:
  ComputeGaussianGridLocal(class LAMMPS *, int, char **);
  ~ComputeGaussianGridLocal() override;
  void init() override;
  void compute_local() override;
  double memory_usage() override;

 protected:
  int ncoeff;
  double **cutsq;
  double rcutfac;     // global cut-off scale
  double *radelem;    // cut-off radius of each atom type
  double *sigmaelem;  // Gaussian width of each atom type
  double *prefacelem; // Gaussian prefactor of each atom type
  double *argfacelem; // Gaussian argument factor of each atom type
  int *map;    // map types to [0,nelements)
  int nelements;
  double cutmax;
};

}    // namespace LAMMPS_NS

#endif
#endif
