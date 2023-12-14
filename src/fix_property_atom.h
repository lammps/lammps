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

#ifdef FIX_CLASS
// clang-format off
FixStyle(property/atom,FixPropertyAtom);
// clang-format on
#else

#ifndef LMP_FIX_PROPERTY_ATOM_H
#define LMP_FIX_PROPERTY_ATOM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtom : public Fix {
 public:
  FixPropertyAtom(class LAMMPS *, int, char **);
  void post_constructor() override;
  ~FixPropertyAtom() override;
  int setmask() override;
  void init() override;

  enum { MOLECULE, CHARGE, RMASS, TEMPERATURE, HEATFLOW, IVEC, DVEC, IARRAY, DARRAY };

  void read_data_section(char *, int, char *, tagint) override;
  bigint read_data_skip_lines(char *) override;
  void write_data_section_size(int, int &, int &) override;
  void write_data_section_pack(int, double **) override;
  void write_data_section_keyword(int, FILE *) override;
  void write_data_section(int, FILE *, int, double **, int) override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_border(int, int *, double *) override;
  int unpack_border(int, int, double *) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;
  double memory_usage() override;

 protected:
  int nvalue, border;
  int molecule_flag, q_flag, rmass_flag;    // flags for specific fields
  int temperature_flag, heatflow_flag;
  int *styles;     // style of each value, see enum
  int *index;      // indices into atom custom data structs
  int *cols;       // columns per value, for arrays
  char *astyle;    // atom style at instantiation

  int values_peratom;    // # of values per atom, including multiple for arrays
  int nmax_old;          // length of peratom arrays the last time they grew
};

}    // namespace LAMMPS_NS

#endif
#endif
