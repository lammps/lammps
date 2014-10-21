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

#ifdef FIX_CLASS

FixStyle(nphug,FixNPHug)

#else

#ifndef LMP_FIX_NPHUG_H
#define LMP_FIX_NPHUG_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNPHug : public FixNH {
 public:
  FixNPHug(class LAMMPS *, int, char **);
  ~FixNPHug();
  void init();
  void setup(int);
  int modify_param(int, char **);
  int pack_restart_data(double *); // pack restart data
  void restart(char *);

 private:
  class Compute *pe;               // PE compute pointer

  void compute_temp_target();
  double compute_vector(int);
  double compute_etotal();
  double compute_vol();
  double compute_hugoniot();
  double compute_us();
  double compute_up();

  char *id_pe;
  int peflag;
  int v0_set,p0_set,e0_set;
  double v0,p0,e0,rho0;
  int idir;
  int uniaxial;

  int size_restart_global();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pstart and Pstop must have the same value

Self-explanatory.

E: Specified target stress must be uniaxial or hydrostatic

Self-explanatory.

E: For triclinic deformation, specified target stress must be hydrostatic

Triclinic pressure control is allowed using the tri keyword, but
non-hydrostatic pressure control can not be used in this case.

E: Temperature control must be used with fix nphug

The temp keyword must be provided.

E: Pressure control must be used with fix nphug

A pressure control keyword (iso, aniso, tri, x, y, or z) must be
provided.

E: Potential energy ID for fix nvt/nph/npt does not exist

A compute for potential energy must be defined.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
