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
DihedralStyle(charmmfsw,DihedralCharmmfsw);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_CHARMMFSW_H
#define LMP_DIHEDRAL_CHARMMFSW_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralCharmmfsw : public Dihedral {
 public:
  DihedralCharmmfsw(class LAMMPS *);
  ~DihedralCharmmfsw() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  int implicit, weightflag, dihedflag;
  double cut_lj_inner14, cut_lj14, cut_coul14;
  double evdwl14_12, evdwl14_6, cut_coulinv14;
  double cut_lj_inner3inv, cut_lj_inner6inv, cut_lj3inv, cut_lj6inv;

  double *k, *weight, *cos_shift, *sin_shift;
  int *multiplicity, *shift;
  double **lj14_1, **lj14_2, **lj14_3, **lj14_4;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
