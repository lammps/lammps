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

FixStyle(CAC/FLUXCHECK,CACFixFluxCheck)

#else

#ifndef LMP_FIX_CAC_FLUXCHECK_H
#define LMP_FIX_CAC_FLUXCHECK_H

#include "fix.h"

namespace LAMMPS_NS {

class CACFixFluxCheck : public Fix {
 public:
  CACFixFluxCheck(class LAMMPS *, int, char **);
  virtual ~CACFixFluxCheck();
  int setmask();
  virtual void init();
  virtual void pre_force(int);
  virtual void setup_pre_force(int);

 protected:
  int flux_compute;
  int ndump;                   // # of Dumps defined
  bigint next_dump_any;        // next timestep for any Dump
  bigint *next_dump;           // next timestep to do each Dump
  bigint *last_dump;           // last timestep each snapshot was output
  class Dump **dump;           // list of defined Dumps
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region for fix oneway does not exist

Self-explanatory.

*/
