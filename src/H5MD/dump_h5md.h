/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   Contributing author: Pierre de Buyl (KU Leuven)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(h5md,DumpH5MD);
// clang-format on
#else

#ifndef LMP_DUMP_H5MD_H
#define LMP_DUMP_H5MD_H

#include "ch5md.h"
#include "dump.h"

namespace LAMMPS_NS {

class DumpH5MD : public Dump {
 public:
  DumpH5MD(class LAMMPS *, int, char **);
  ~DumpH5MD() override;

 private:
  int natoms, ntotal;
  int unwrap_flag;    // 1 if atom coords are unwrapped, 0 if no
  h5md_file datafile;
  h5md_particles_group particles_data;
  char *author_name;
  DumpH5MD *other_dump;

  bool do_box;
  bool create_group;

  // data arrays and intervals
  int every_dump;
  double *dump_position;
  int every_position;
  int *dump_image;
  int every_image;
  double *dump_velocity;
  int every_velocity;
  double *dump_force;
  int every_force;
  int *dump_species;
  int every_species;
  int *dump_charge;
  int every_charge;

  void init_style() override;
  int modify_param(int, char **) override;
  void openfile() override;
  void write_header(bigint) override{};
  void pack(tagint *) override;
  void write_data(int, double *) override;

  void write_frame();
  void write_fixed_frame();
};

}    // namespace LAMMPS_NS

#endif
#endif
