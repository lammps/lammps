/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
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
  virtual void force_zeta(Param *, double, double, double &, double &, int, double &);

};

}

#endif
#endif
