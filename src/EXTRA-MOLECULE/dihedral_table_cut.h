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
DihedralStyle(table/cut,DihedralTableCut);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_TABLE_CUT_H
#define LMP_DIHEDRAL_TABLE_CUT_H

#include "dihedral_table.h"

namespace LAMMPS_NS {

class DihedralTableCut : public DihedralTable {
 public:
  DihedralTableCut(class LAMMPS *);
  ~DihedralTableCut() override;
  void compute(int, int) override;
  void coeff(int, char **) override;

 protected:
  double *aat_k, *aat_theta0_1, *aat_theta0_2;

  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
