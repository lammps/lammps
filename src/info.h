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

#ifdef COMMAND_CLASS

CommandStyle(info,Info)

#else

#ifndef LMP_INFO_H
#define LMP_INFO_H

#include "pointers.h"

namespace LAMMPS_NS {

class Info : protected Pointers {
 public:
  Info(class LAMMPS *lmp) : Pointers(lmp) {};
  void command(int, char **);

  bool is_active(const char *, const char *);
  bool is_defined(const char *, const char *);
  bool is_available(const char *, const char *);

  static bool has_gzip_support();
  static bool has_png_support();
  static bool has_jpeg_support();
  static bool has_ffmpeg_support();
  static bool has_exceptions();
  static bool has_package(const char * package_name);

  static char *get_os_info();
  static char *get_compiler_info();

  char **get_variable_names(int &num);

private:
  void available_styles(FILE * out, int flags);

  void atom_styles(FILE * out);
  void integrate_styles(FILE * out);
  void minimize_styles(FILE * out);
  void pair_styles(FILE * out);
  void bond_styles(FILE * out);
  void angle_styles(FILE * out);
  void dihedral_styles(FILE * out);
  void improper_styles(FILE * out);
  void kspace_styles(FILE * out);
  void fix_styles(FILE * out);
  void compute_styles(FILE * out);
  void region_styles(FILE * out);
  void dump_styles(FILE * out);
  void command_styles(FILE * out);
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Ignoring unknown or incorrect info command flag

Self-explanatory.  An unknown argument was given to the info command.
Compare your input with the documentation.

E: Unknown name for info package category

Self-explanatory.

E: Unknown name for info newton category

Self-explanatory.

E: Unknown name for info pair category

Self-explanatory.

E: Unknown category for info is_active()

Self-explanatory.

E: Unknown category for info is_available()

Self-explanatory.

E: Unknown category for info is_defined()

Self-explanatory.

*/
