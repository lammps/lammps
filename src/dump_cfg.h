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

DumpStyle(cfg,DumpCFG)

#else

#ifndef LMP_DUMP_CFG_H
#define LMP_DUMP_CFG_H

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpCFG : public DumpCustom {
 public:
  int multifile_override;          // used by write_dump command

  DumpCFG(class LAMMPS *, int, char **);
  virtual ~DumpCFG();

 private:
  char **auxname;            // name strings of auxiliary properties
  int unwrapflag;            // 1 if unwrapped coordinates are requested

  void init_style();
  void write_header(bigint);
  int convert_string(int, double *);
  void write_data(int, double *);

  typedef void (DumpCFG::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;             // ptr to write data functions
  void write_string(int, double *);
  void write_lines(int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Dump cfg arguments must start with 'mass type xs ys zs' or 'mass type xsu ysu zsu'

This is a requirement of the CFG output format.  See the dump cfg doc
page for more details.

E: Dump cfg arguments can not mix xs|ys|zs with xsu|ysu|zsu

Self-explanatory.

E: Invalid keyword in dump cfg command

Self-explanatory.

E: Dump cfg requires one snapshot per file

Use the wildcard "*" character in the filename.

*/
