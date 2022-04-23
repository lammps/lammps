/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(dcd,DumpDCD);
// clang-format on
#else

#ifndef LMP_DUMP_DCD_H
#define LMP_DUMP_DCD_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpDCD : public Dump {
 public:
  DumpDCD(LAMMPS *, int, char **);
  ~DumpDCD() override;

 private:
  int natoms, ntotal;
  int headerflag, nevery_save, nframes;

  float *coords, *xf, *yf, *zf;
  int unwrap_flag;    // 1 if atom coords are unwrapped, 0 if no

  void init_style() override;
  void openfile() override;
  void write_header(bigint) override;
  void pack(tagint *) override;
  void write_data(int, double *) override;
  int modify_param(int, char **) override;
  double memory_usage() override;

  void write_frame();
  void write_dcd_header(const char *);
};

}    // namespace LAMMPS_NS

#endif
#endif
