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

DumpStyle(custom/time,DumpCustom_Time)

#else

#ifndef LMP_DUMP_CUSTOM_TIME_H
#define LMP_DUMP_CUSTOM_TIME_H

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpCustom_Time : public DumpCustom {
 public:
  DumpCustom_Time(class LAMMPS *, int, char **);
  virtual ~DumpCustom_Time();

 protected:
  double time_every;           //time frequency of a dump
  double next_time;            //next time for a dump
  double tol;        
  int write_time;              // -1-write in the next timestep, 0-write now, 1-check if it's time for a write
  int nvarcompute;             //number of computes used by variables in dump
  class Compute ** varcompute; //list of ptrs to the Compute objects used in variables in dump

  // private methods

  virtual void write_header(bigint);
  virtual void write_data(int, double *);
  int count();
  virtual int modify_param(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

The same as in dump_custom.h

*/
