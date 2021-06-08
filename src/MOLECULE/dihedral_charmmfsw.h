/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
  virtual ~DihedralCharmmfsw();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  virtual void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

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

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect multiplicity arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect weight arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Dihedral style charmmfsw must be set to same r-RESPA level as 'pair'

UNDOCUMENTED

E: Dihedral style charmmfsw must be set to same r-RESPA level as 'outer'

UNDOCUMENTED

E: Must use 'special_bonds charmm' with dihedral style charmm for use with CHARMM pair styles

UNDOCUMENTED

E: Dihedral charmmfsw is incompatible with Pair style

Dihedral style charmmfsw must be used with a pair style charmm
in order for the 1-4 epsilon/sigma parameters to be defined.

*/
