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

PairStyle(tersoff/zbl/omp,PairTersoffZBLOMP)

#else

#ifndef LMP_PAIR_TERSOFF_ZBL_OMP_H
#define LMP_PAIR_TERSOFF_ZBL_OMP_H

#include "pair_tersoff_omp.h"

namespace LAMMPS_NS {

class PairTersoffZBLOMP : public PairTersoffOMP {
 public:
  PairTersoffZBLOMP(class LAMMPS *);
  ~PairTersoffZBLOMP() {}

 private:
  double global_a_0;		// Bohr radius for Coulomb repulsion
  double global_epsilon_0;	// permittivity of vacuum for Coulomb repulsion
  double global_e;		// proton charge (negative of electron charge)

  void read_file(char *);
  void repulsive(const Param &, const double &, double &, int, double &);

  double ters_fa(const double &, const Param &) const;
  double ters_fa_d(const double &, const Param &) const;
	
/* ----------------------------------------------------------------------
   Fermi-like smoothing function
------------------------------------------------------------------------- */

  double F_fermi(const double &r, const Param &param) const {
    return 1.0 / (1.0 + exp(-param.ZBLexpscale*(r-param.ZBLcut)));
  };

/* ----------------------------------------------------------------------
   Fermi-like smoothing function derivative with respect to r
------------------------------------------------------------------------- */

  double F_fermi_d(const double &r, const Param &param) const {
    double tmp = 1.0 + exp(-param.ZBLexpscale*(r-param.ZBLcut));

    return param.ZBLexpscale*exp(-param.ZBLexpscale*(r-param.ZBLcut)) 
      / (tmp*tmp);
  };
};

}

#endif
#endif
