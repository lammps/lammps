/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
  virtual ~DumpH5MD();

 private:
  int natoms, ntotal;
  int unwrap_flag;    // 1 if atom coords are unwrapped, 0 if no
  h5md_file datafile;
  int datafile_from_dump;
  h5md_particles_group particles_data;
  char *author_name;

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

  void init_style();
  int modify_param(int, char **);
  void openfile();
  void write_header(bigint);
  void pack(tagint *);
  void write_data(int, double *);

  void write_frame();
  void write_fixed_frame();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid number of arguments in dump h5md

Make sure that each data item (position, etc.) is followed by a dump
interval.

E: Dump h5md requires sorting by atom ID

Use the dump_modify sort command to enable this.

E: Cannot use variable every setting for dump xtc

The format of this file requires snapshots at regular intervals.

E: Cannot change dump_modify every for dump xtc

The frequency of writing dump xtc snapshots cannot be changed.

E: Cannot set file_from in dump h5md after box or create_group

The file_from option modifies the box and create_group options and
they must appear after file_from if used.

*/
