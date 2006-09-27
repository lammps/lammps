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

#ifndef DUMP_DCD_H
#define DUMP_DCD_H

#include "stdio.h"
#include "dump.h"
#include "inttypes.h"

class DumpDCD : public Dump {
 public:
  DumpDCD(int, char**);
  ~DumpDCD();
  void init();
  int memory_usage();

 private:
  int natoms,ntotal,headerflag,nevery_save,nframes;
  float *coords,*xf,*yf,*zf;

  void openfile();
  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);

  void write_frame();
  void write_dcd_header(char *);
  void fwrite_int32(FILE *, uint32_t);
};

#endif
