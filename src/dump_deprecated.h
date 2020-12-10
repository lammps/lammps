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

#ifdef DUMP_CLASS

// list all deprecated and removed dump styles here

DumpStyle(DEPRECATED,DumpDeprecated)

#else

#ifndef LMP_DUMP_DEPRECATED_H
#define LMP_DUMP_DEPRECATED_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpDeprecated : public Dump {
 public:
  DumpDeprecated(class LAMMPS *, int, char **);
  ~DumpDeprecated() {}
  virtual void init_style() {}
  virtual void write_header(bigint) {}
  virtual void pack(tagint *) {}
  virtual void write_data(int, double *) {}
 };

}

#endif
#endif

/* ERROR/WARNING messages:

E: This dump style has been removed from LAMMPS

UNDOCUMENTED

*/
