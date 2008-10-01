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

#ifndef PAIR_TERSOFF_ZBL_H
#define PAIR_TERSOFF_ZBL_H

#include "pair_tersoff.h"

namespace LAMMPS_NS {

class PairTersoffZBL : public PairTersoff {
 public:
  PairTersoffZBL(class LAMMPS *);
  ~PairTersoffZBL() {}

 private:
  double global_a_0;		// Bohr radius for Coulomb repulsion
  double global_epsilon_0;	// permittivity of vacuum for Coulomb repulsion
  double global_e;		// proton charge (negative of electron charge)

  void read_file(char *);
  void repulsive(Param *, double, double &, int, double &);

  double ters_fa(double, Param *);
  double ters_fa_d(double, Param *);
	
  double F_fermi(double, Param *);
  double F_fermi_d(double, Param *);
};

}

#endif
