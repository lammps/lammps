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
  virtual ~PairTersoffZBLOMP() {}

 protected:
  double global_a_0;		// Bohr radius for Coulomb repulsion
  double global_epsilon_0;	// permittivity of vacuum for Coulomb repulsion
  double global_e;		// proton charge (negative of electron charge)

  virtual void read_file(char *);
  virtual void repulsive(Param *, double, double &, int, double &);

  virtual double ters_fa(double, Param *);
  virtual double ters_fa_d(double, Param *);
	
  virtual double F_fermi(double, Param *);
  virtual double F_fermi_d(double, Param *);
};

}

#endif
#endif
