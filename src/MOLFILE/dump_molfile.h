/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(molfile,DumpMolfile);
// clang-format on
#else

#ifndef LMP_DUMP_MOLFILE_H
#define LMP_DUMP_MOLFILE_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpMolfile : public Dump {
 public:
  DumpMolfile(LAMMPS *, int, char **);
  ~DumpMolfile() override;
  void write() override;

 protected:
  class MolfileInterface *mf;    //< handles low-level I/O
  // per-atom data
  float *coords, *vels, *masses, *charges, *radiuses;
  int *types, *molids;
  char **typenames;

  int natoms, me, ntotal, ntypes;
  int need_structure;
  int unwrap_flag;      // 1 if writing unwrapped atom coords, 0 if not
  int velocity_flag;    // 1 if writing velocities, 0 if not
  int topology_flag;    // 1 if writing topology data, 0 if not
  float cell[6];        // cell parameters: A, B, C, alpha, beta, gamma

  void init_style() override;
  int modify_param(int, char **) override;
  void write_header(bigint) override{};
  void pack(tagint *) override;
  void write_data(int, double *) override;
  double memory_usage() override;
  void openfile() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
