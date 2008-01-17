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

#ifndef DUMP_XTC_H
#define DUMP_XTC_H

#include "dump.h"

#ifdef LAMMPS_XDR
#include "xdr_compat.h"
#else
#include "rpc/rpc.h"
#include "rpc/xdr.h"
#endif

namespace LAMMPS_NS {

class DumpXTC : public Dump {
 public:
  DumpXTC(class LAMMPS *, int, char**);
  ~DumpXTC();
  void init();
  double memory_usage();
	
 private:
  int natoms,ntotal;
  int unwrap_flag;            // 1 if atom coords are unwrapped, 0 if no
  float precision;            // user-adjustable precision setting
  float *coords;
  double sfactor;
  XDR xd;

  int modify_param(int, char **);
  void openfile();
  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);

  void write_frame();
};

}

#endif
