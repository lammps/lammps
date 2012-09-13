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

PairStyle(born/coul/msm,PairBornCoulMSM)

#else

#ifndef LMP_PAIR_BORN_COUL_MSM_H
#define LMP_PAIR_BORN_COUL_MSM_H

#include "pair_born_coul_long.h"

namespace LAMMPS_NS {

class PairBornCoulMSM : public PairBornCoulLong {
 public:
  PairBornCoulMSM(class LAMMPS *);
  virtual ~PairBornCoulMSM(){};
  virtual void compute(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

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

E: Pair style born/coul/long requires atom attribute q

An atom style that defines this attribute must be used.

E: Pair style is incompatible with KSpace style

If a pair style with a long-range Coulombic component is selected,
then a kspace style must also be used.

*/
