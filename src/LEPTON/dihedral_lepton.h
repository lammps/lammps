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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(lepton,DihedralLepton);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_LEPTON_H
#define LMP_DIHEDRAL_LEPTON_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralLepton : public Dihedral {
 public:
  DihedralLepton(class LAMMPS *);
  ~DihedralLepton() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  std::vector<std::string> expressions;
  int *type2expression;

  virtual void allocate();
  double get_phi(double const *, double const *, double const *, double const *,
                 class Domain *domain, double *, double *, double *, double *, double *) const;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval();
};
}    // namespace LAMMPS_NS
#endif
#endif
