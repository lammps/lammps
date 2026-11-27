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
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#ifndef LMP_LDD_POTENTIAL_H
#define LMP_LDD_POTENTIAL_H

#include "pointers.h"

namespace LAMMPS_NS {

class LddPotential : protected Pointers {
 public:
  int allocated;      // 0 or 1, tracks whether the potential type has been allocated
  int n_coeffs;       // length of coeffs
  double *coeffs;     // coefficients involved in defining u_LD or u_SG (e.g. a b c in ax^2 + bx +c)

  struct t_table {
    int n_pts;     // number of table entries
    double dr;     // domain spacing between table entries
    double *r;     // the x column list of entries from the table (domain) (usually read in)
    double *u;     // the y column list of entries from the table (range) (usually read in)
    double *u2;    // the second derivative of the y entries at the tabulated knots (inferred)
    double *
        f;    // -dU/dr , either read from the third column of the table or inferred from the second (see table/spline vs. table/gradspline
    double *f2;    // The second derivatives calculated from f (inferred)
  };    // structure to hold tabulated table info

  t_table potl_table;    // field for the table info

  LddPotential(class LAMMPS *);
  ~LddPotential() override;

  virtual void setup_potl(int, int, char **) = 0;    // fnc to define U_x
  virtual double u(double) = 0;                      // fn that returns value of U_x(rho)
  virtual double f(double) = 0;                      // fn that returns f_x(rho)

  void read_table_file(char *, bool);
  int get_table_index(double);
  double calc_A_table(double, int);
};

}    // namespace LAMMPS_NS

#endif
