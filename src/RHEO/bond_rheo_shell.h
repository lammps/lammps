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

#ifdef BOND_CLASS
// clang-format off
BondStyle(rheo/shell,BondRHEOShell);
// clang-format on
#else

#ifndef LMP_BOND_RHEO_SHELL_H
#define LMP_BOND_RHEO_SHELL_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondRHEOShell : public BondBPM {
 public:
  BondRHEOShell(class LAMMPS *);
  ~BondRHEOShell() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *k, *ecrit, *gamma;
  double tform, rmax, rsurf;

  int *dbond, *nbond;
  int index_nb, nmax_store;
  char *id_fix;

  class ComputeRHEOSurface *compute_surface;

  void process_ineligibility(int, int);
  void allocate();
  void store_data();
  double store_bond(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
