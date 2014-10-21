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

PairStyle(buck/coul/msm,PairBuckCoulMSM)

#else

#ifndef LMP_PAIR_BUCK_COUL_MSM_H
#define LMP_PAIR_BUCK_COUL_MSM_H

#include "pair_buck_coul_long.h"

namespace LAMMPS_NS {

class PairBuckCoulMSM : public PairBuckCoulLong {
 public:
  PairBuckCoulMSM(class LAMMPS *);
  virtual ~PairBuckCoulMSM();
  virtual void compute(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);
 protected:
  int nmax;
  double **ftmp;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Pair style buck/coul/long requires atom attribute q

The atom style defined does not have these attributes.

E: Pair style is incompatible with KSpace style

If a pair style with a long-range Coulombic component is selected,
then a kspace style must also be used.

E: Must use 'kspace_modify pressure/scalar no' to obtain per-atom virial with kspace_style MSM

The kspace scalar pressure option cannot be used to obtain per-atom virial.

*/
