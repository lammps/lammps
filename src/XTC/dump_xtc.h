/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef DUMP_XTC_H
#define DUMP_XTC_H

#include "dump.h"
#include "rpc/rpc.h"
#include "rpc/xdr.h"

class DumpXTC : public Dump {
 public:
  DumpXTC(int, char**);
  ~DumpXTC();
  void init();
  int memory_usage();
	
 private:
  int natoms,ntotal;
  float precision;
  float *coords;
  double sfactor;
  XDR xd;

  void openfile();
  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);

  void write_frame();
};

#endif
