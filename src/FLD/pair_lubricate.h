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

PairStyle(lubricate,PairLubricate)

#else

#ifndef LMP_PAIR_LUBRICATE_H
#define LMP_PAIR_LUBRICATE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLubricate : public Pair {
 public:
  PairLubricate(class LAMMPS *);
  virtual ~PairLubricate();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  virtual void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  int pre_adapt(char *, int, int, int, int);
  void adapt(int, int, int, int, int, double);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 protected:
  double mu,cut_inner_global,cut_global;
  double rad;
  int flaglog,flagfld,shearing;
  int flagdeform, flagwall;
  double vol_P;
  class FixWall *wallfix;
  int flagVF, flagHI;

  double Ef[3][3];
  double R0,RT0,RS0;
  double **cut_inner,**cut;

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

E: Pair lubricate requires atom style sphere

Self-explanatory.

E: Pair lubricate requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Pair lubricate requires monodisperse particles

All particles must be the same finite size.

E: Using pair lubricate with inconsistent fix deform remap option

Must use remap v option with fix deform with this pair style.

E: Cannot use multiple fix wall commands with pair lubricate

Self-explanatory.

*/
