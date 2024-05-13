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

#ifdef FIX_CLASS
// clang-format off
FixStyle(indent,FixIndent);
// clang-format on
#else

#ifndef LMP_FIX_INDENT_H
#define LMP_FIX_INDENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixIndent : public Fix {
 public:
  FixIndent(class LAMMPS *, int, char **);
  ~FixIndent() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  int istyle, scaleflag, side;
  double k, k3;
  char *xstr, *ystr, *zstr, *rstr, *pstr;
  int xvar, yvar, zvar, rvar, pvar;
  double xvalue, yvalue, zvalue, rvalue, pvalue;
  int indenter_flag, planeside;
  double indenter[4], indenter_all[4];
  int cdim, varflag;
  int ilevel_respa;

  char *rlostr, *rhistr, *lostr, *histr;
  int rlovar, rhivar, lovar, hivar;
  double rlovalue, rhivalue, lovalue, hivalue;

  // methods for argument parsing

  int geometry(int, char **);
  void options(int, char **);

  // methods for conical indenter

  bool PointInsideCone(int, double *, double, double, double, double, double *);
  void DistanceExteriorPoint(int, double *, double, double, double, double,
                             double &, double &, double &);
  void DistanceInteriorPoint(int, double *, double, double, double, double,
                             double &, double &, double &);
  void point_on_line_segment(double *, double *, double *, double *);
  double closest(double *, double *, double *, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
