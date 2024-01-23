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
BondStyle(oxdna/fene,BondOxdnaFene);
// clang-format on
#else

#ifndef LMP_BOND_OXDNA_FENE_H
#define LMP_BOND_OXDNA_FENE_H

#include "bond.h"

namespace LAMMPS_NS {

class BondOxdnaFene : public Bond {
 public:
  BondOxdnaFene(class LAMMPS *lmp) : Bond(lmp) {}
  ~BondOxdnaFene() override;
  virtual void compute_interaction_sites(double *, double *, double *, double *) const;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *k, *Delta, *r0;                       // FENE
  double **nx_xtrct, **ny_xtrct, **nz_xtrct;    // per-atom arrays for local unit vectors

  void allocate();
  void ev_tally_xyz(int, int, int, int, double, double, double, double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
