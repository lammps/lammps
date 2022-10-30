/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(cfg,DumpCFG);
// clang-format on
#else

#ifndef LMP_DUMP_CFG_H
#define LMP_DUMP_CFG_H

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpCFG : public DumpCustom {
 public:
  int multifile_override;    // used by write_dump command

  DumpCFG(class LAMMPS *, int, char **);
  ~DumpCFG() override;

 protected:
  char **auxname;    // name strings of auxiliary properties
  int unwrapflag;    // 1 if unwrapped coordinates are requested

  void init_style() override;
  void write_header(bigint) override;
  int convert_string(int, double *) override;
  void write_data(int, double *) override;

  typedef void (DumpCFG::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_string(int, double *);
  void write_lines(int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
