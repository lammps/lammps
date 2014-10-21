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

PairStyle(brownian,PairBrownian)

#else

#ifndef LMP_PAIR_BROWNIAN_H
#define LMP_PAIR_BROWNIAN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBrownian : public Pair {
 public:
  PairBrownian(class LAMMPS *);
  virtual ~PairBrownian();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  double cut_inner_global,cut_global;
  double t_target,mu;
  int flaglog,flagfld;
  int flagHI, flagVF;
  int flagdeform, flagwall;
  double vol_P;
  double rad;
  class FixWall *wallfix;

  int seed;
  double **cut_inner,**cut;
  double R0,RT0;

  class RanMars *random;

  void set_3_orthogonal_vectors(double*,double*,double*);
  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Cannot include log terms without 1/r terms; setting flagHI to 1

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair brownian requires atom style sphere

Self-explanatory.

W: Pair brownian needs newton pair on for momentum conservation

Self-explanatory.

E: Pair brownian requires extended particles

One of the particles has radius 0.0.

E: Pair brownian requires monodisperse particles

All particles must be the same finite size.

E: Cannot use multiple fix wall commands with pair brownian

Self-explanatory.

*/
