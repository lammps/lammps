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

#ifndef LMP_VARIABLE_KOKKOS_H
#define LMP_VARIABLE_KOKKOS_H

#include "variable.h"

namespace LAMMPS_NS {

class VariableKokkos : public Variable {
 public:
  VariableKokkos(class LAMMPS *lmp) : Variable(lmp) {}

  void compute_atom(int, int, double *, int, int) override;

 protected:
  void atom_vector(char *, Tree **, Tree **, int &) override;
  int group_function(char *, char *, Tree **, Tree **, int &, double *, int &, int) override;
  int special_function(const std::string &, char *, Tree **, Tree **, int &, double *, int &, int,
                       char *, int &, char *&) override;
  void peratom2global(int, char *, double *, int, tagint, Tree **, Tree **, int &, double *,
                      int &) override;
  void custom2global(int *, double *, int, tagint, Tree **, Tree **, int &, double *,
                     int &) override;
  void sync_peratom(const char *) override;

 private:
  void sync_host(uint64_t);
};

}    // namespace LAMMPS_NS

#endif
