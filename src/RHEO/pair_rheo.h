/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(rheo,PairRHEO)

#else

#ifndef LMP_PAIR_RHEO_H
#define LMP_PAIR_RHEO_H

#include "pair.h"

namespace LAMMPS_NS {

class PairRHEO : public Pair {
 public:
  PairRHEO(class LAMMPS *);
  ~PairRHEO() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void setup() override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  int laplacian_order;          // From fix RHEO
  double h, csq*, rho0*;        // From fix RHEO

  double hsq, hinv, rho0, av, rho_damp;

  int artificial_visc_flag;
  int rho_damp_flag;
  int thermal_flag;

  void allocate();

  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOGrad *compute_grad;
  class ComputeRHEOSolidInterpolation *compute_sinterpolation;
  class FixRHEO *fix_rheo;
};

}

#endif
#endif
