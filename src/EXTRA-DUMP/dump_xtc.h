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
DumpStyle(xtc,DumpXTC);
// clang-format on
#else

#ifndef LMP_DUMP_XTC_H
#define LMP_DUMP_XTC_H

#include "dump.h"

struct XDR;
namespace LAMMPS_NS {

class DumpXTC : public Dump {
 public:
  DumpXTC(class LAMMPS *, int, char **);
  ~DumpXTC() override;

 private:
  int natoms, ntotal;
  int nevery_save;
  int unwrap_flag;    // 1 if atom coords are unwrapped, 0 if no
  float precision;    // user-adjustable precision setting
  float *coords;
  double sfactor, tfactor;    // scaling factors for positions and time unit
  XDR *xd;

  void init_style() override;
  int modify_param(int, char **) override;
  void openfile() override;
  void write_header(bigint) override;
  void pack(tagint *) override;
  void write_data(int, double *) override;
  double memory_usage() override;

  void write_frame();
};

}    // namespace LAMMPS_NS

#endif
#endif
