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
DihedralStyle(class2,DihedralClass2);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_CLASS2_H
#define LMP_DIHEDRAL_CLASS2_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralClass2 : public Dihedral {
 public:
  DihedralClass2(class LAMMPS *);
  ~DihedralClass2() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  double *k1, *k2, *k3;
  double *phi1, *phi2, *phi3;
  double *mbt_f1, *mbt_f2, *mbt_f3, *mbt_r0;
  double *ebt_f1_1, *ebt_f2_1, *ebt_f3_1, *ebt_r0_1;
  double *ebt_f1_2, *ebt_f2_2, *ebt_f3_2, *ebt_r0_2;
  double *at_f1_1, *at_f2_1, *at_f3_1, *at_theta0_1;
  double *at_f1_2, *at_f2_2, *at_f3_2, *at_theta0_2;
  double *aat_k, *aat_theta0_1, *aat_theta0_2;
  double *bb13t_k, *bb13t_r10, *bb13t_r30;
  int *setflag_d, *setflag_mbt, *setflag_ebt;
  int *setflag_at, *setflag_aat, *setflag_bb13t;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
