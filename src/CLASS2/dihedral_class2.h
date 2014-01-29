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

#ifdef DIHEDRAL_CLASS

DihedralStyle(class2,DihedralClass2)

#else

#ifndef LMP_DIHEDRAL_CLASS2_H
#define LMP_DIHEDRAL_CLASS2_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralClass2 : public Dihedral {
 public:
  DihedralClass2(class LAMMPS *);
  virtual ~DihedralClass2();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  double *k1,*k2,*k3;
  double *phi1,*phi2,*phi3;
  double *mbt_f1,*mbt_f2,*mbt_f3,*mbt_r0;
  double *ebt_f1_1,*ebt_f2_1,*ebt_f3_1,*ebt_r0_1;
  double *ebt_f1_2,*ebt_f2_2,*ebt_f3_2,*ebt_r0_2;
  double *at_f1_1,*at_f2_1,*at_f3_1,*at_theta0_1;
  double *at_f1_2,*at_f2_2,*at_f3_2,*at_theta0_2;
  double *aat_k,*aat_theta0_1,*aat_theta0_2;
  double *bb13t_k,*bb13t_r10,*bb13t_r30;
  int *setflag_d,*setflag_mbt,*setflag_ebt;
  int *setflag_at,*setflag_aat,*setflag_bb13t;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld    

UNDOCUMENTED

E: Invalid coeffs for this dihedral style

Cannot set class 2 coeffs in data file for this dihedral style.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

U: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

*/
