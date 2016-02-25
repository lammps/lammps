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

DumpStyle(cfg/gz,DumpCFGGZ)

#else

#ifndef LMP_DUMP_CFG_GZ_H
#define LMP_DUMP_CFG_GZ_H

#include "dump_cfg.h"
#include <zlib.h>

namespace LAMMPS_NS {

class DumpCFGGZ : public DumpCFG {
 public:
  DumpCFGGZ(class LAMMPS *, int, char **);
  virtual ~DumpCFGGZ();

 protected:
  gzFile gzFp;  // file pointer for the compressed output stream

  virtual void openfile();
  virtual void write_header(bigint);
  virtual void write_data(int, double *);
  virtual void write();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Dump cfg/gz only writes compressed files

The dump cfg/gz output file name must have a .gz suffix.

E: Cannot open dump file

Self-explanatory.

*/
