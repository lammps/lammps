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

#ifndef DUMP_H
#define DUMP_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Dump : protected Pointers {
 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump
  int igroup,groupbit;       // group that Dump is performed on
  int me,nprocs;             // proc info

  char *filename;            // user-specified file
  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all, 1 = one file/proc

  int header_flag;           // 0 = item, 2 = xyz
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int sort_flag;             // 1 if write in sorted order, 0 if not
  int append_flag;           // 1 if open file in append mode, 0 if not
  int singlefile_opened;     // 1 = one big file, already opened, else 0

  char *format_default;      // default format string
  char *format_user;         // format string set by user
  char *format;              // format string for the file write
  double *buf;               // memory for atom quantities
  int maxbuf;                // size of buf
  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one atom

  Dump(class LAMMPS *, int, char **);
  virtual ~Dump();
  virtual void init() {}
  void write();
  void modify_params(int, char **);
  virtual double memory_usage();

 protected:
  double boxxlo,boxxhi;      // local copies of domain values
  double boxylo,boxyhi;      // lo/hi are bounding box for triclinic
  double boxzlo,boxzhi;
  double boxxy,boxxz,boxyz;

  virtual void openfile();
  virtual int modify_param(int, char **) {return 0;}

  virtual void write_header(int) = 0;
  virtual int count() = 0;
  virtual int pack() = 0;
  virtual void write_data(int, double *) = 0;
};

}

#endif
