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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mesocnt, PairMesoCNT);
// clang-format on
#else

#ifndef LMP_PAIR_MESOCNT_H
#define LMP_PAIR_MESOCNT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMesoCNT : public Pair {
 public:
  PairMesoCNT(class LAMMPS *);
  ~PairMesoCNT();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  int uinf_points, gamma_points, phi_points, usemi_points;
  int nlocal_size, reduced_neigh_size;
  int *reduced_nlist, *numchainlist;
  int **reduced_neighlist, **nchainlist, **endlist;
  int ***chainlist;

  double ang, ang_inv, eunit, funit;
  double delta1, delta2;
  double r, rsq, d, rc, rcsq, rc0, cutoff, cutoffsq;
  double r_ang, rsq_ang, d_ang, rc_ang, rcsq_ang, cutoff_ang, cutoffsq_ang;
  double sig, sig_ang, comega, ctheta;
  double hstart_uinf, hstart_gamma, hstart_phi, psistart_phi, hstart_usemi, xistart_usemi;
  double delh_uinf, delh_gamma, delh_phi, delpsi_phi, delh_usemi, delxi_usemi;

  double p1[3], p2[3], p[3], m[3];
  double *param, *w, *wnode;
  double **dq_w;
  double ***q1_dq_w, ***q2_dq_w;
  double *uinf_data, *gamma_data, **phi_data, **usemi_data;
  double **uinf_coeff, **gamma_coeff, ****phi_coeff, ****usemi_coeff;
  double **flocal, **fglobal, **basis;

  char *file;

  void allocate();
  void bond_neigh();
  void neigh_common(int, int, int &, int *);
  void chain_split(int *, int, int &, int **, int *, int *);
  void sort(int *, int);
  void read_file();
  void read_data(FILE *, double *, double &, double &, int);
  void read_data(FILE *, double **, double &, double &, double &, double &, int);

  void spline_coeff(double *, double **, double, int);
  void spline_coeff(double **, double ****, double, double, int);

  double spline(double, double, double, double **, int);
  double dspline(double, double, double, double **, int);
  double spline(double, double, double, double, double, double, double ****, int);
  double dxspline(double, double, double, double, double, double, double ****, int);
  double dyspline(double, double, double, double, double, double, double ****, int);

  void geometry(const double *, const double *, const double *, const double *, const double *,
                double *, double *, double *, double **);
  void weight(const double *, const double *, const double *, const double *, double &, double *,
              double *, double *, double *);

  void finf(const double *, double &, double **);
  void fsemi(const double *, double &, double &, double **);

  // inlined functions for efficiency

  inline double heaviside(double x)
  {
    if (x > 0)
      return 1.0;
    else
      return 0.0;
  }

  inline double s(double x)
  {
    return heaviside(-x) + heaviside(x) * heaviside(1 - x) * (1 - x * x * (3 - 2 * x));
  }

  inline double ds(double x) { return 6 * heaviside(x) * heaviside(1 - x) * x * (x - 1); }

  inline double s5(double x)
  {
    double x2 = x * x;
    return heaviside(-x) + heaviside(x) * heaviside(1 - x) * (1 - x2 * x * (6 * x2 - 15 * x + 10));
  }

  inline double ds5(double x)
  {
    double x2 = x * x;
    return -30 * heaviside(x) * heaviside(1 - x) * x2 * (x2 - 2 * x + 1);
  }
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style mesocnt does not support lj units

Self-explanatory. Specify different unit system using the units
command.

E: Pair style mesocnt requires atom IDs

Self-explanatory. Turn on atom IDs using the atom_modify command.

E: Pair style mesocnt requires newton pair on

Self-explanatory. Turn on Newton's third law with the newton command.

E: Cannot open mesocnt file %s

The specified mesocnt potential file cannot be opened. Check that the
path and name are correct.

E: Premature end of file in pair table %s

The specified mesocnt potential file is shorter than specified. Check
if the correct file is being used and the right number of data points
was specified in the pair_style

W: %d of %d lines were incomplete or could not be parsed completely
in pair table %s

A number of lines in the specified mesocnt potential file is incomplete
or in the wrong format. Check the file for errors and missing data.

W: %d spacings in the first column were different from the first spacing
in the pair table %s

The spacings between x coordinates in the first column of the specified
mesocnt potential file vary throughout the file. Use a potential file
with higher precision.

W: %d spacings in second column were different from first
spacing in pair table %s

The spacings between y coordinates in the second column of the specified
mesocnt potential file vary throughout the file. Use a potential file
with higher precision.

*/
