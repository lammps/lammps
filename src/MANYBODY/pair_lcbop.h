/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lcbop,PairLCBOP)

#else

#ifndef LMP_PAIR_LCBOP_H
#define LMP_PAIR_LCBOP_H

#include "pair.h"
#include "my_page.h"
#include "math.h"
#include "math_const.h"

namespace LAMMPS_NS {

class PairLCBOP : public Pair {
 public:
  PairLCBOP(class LAMMPS *);
  virtual ~PairLCBOP();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 protected:
  int **pages;                     // neighbor list pages
  int *map;                        // 0 (C) or -1 (NULL) for each type

  int me;

  double cutLR;                    // LR cutoff

  double cutLRsq;                  // LR cutoff squared
  double cut3rebo;                 // maximum distance for 3rd SR neigh

  int maxlocal;                    // size of numneigh, firstneigh arrays
  int maxpage;                     // # of pages currently allocated
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  MyPage<int> *ipage;              // neighbor list pages
  int *SR_numneigh;                // # of pair neighbors for each atom
  int **SR_firstneigh;             // ptr to 1st neighbor of each atom

  double *N;                       // sum of cutoff fns ( f_C ) with SR neighs
  double *M;                       // sum_j f_C_ij*F(N_j - f_C_ij)

  double
    r_1, r_2, gamma_1, A, B_1, B_2, alpha, beta_1, beta_2,
    d, C_1, C_4, C_6, L, kappa, R_0, R_1,
    r_0, r_1_LR, r_2_LR,
    v_1, v_2, eps_1, eps_2, lambda_1, lambda_2, eps, delta;
  double r_2_sq;

  // splines coefficients
  struct TF_conj_field {
    double
        f_00,
        f_01,
        f_10,
        f_11,
        f_x_00,
        f_x_01,
        f_x_10,
        f_x_11,
        f_y_00,
        f_y_01,
        f_y_10,
        f_y_11;
  } F_conj_field[3][3][2];

  double F_conj_data[4][4][2][3]; // temporary data from file
  double gX[6];        // x coordinates for described points[# of points];
  double gC[5+1][6-1]; // coefficients for each period between described points [degree of polynomial+1][# of points-1]

  void SR_neigh();
  void FSR(int, int);
  void FLR(int, int);

  void FNij( int, int, double, double**, int );
  void FMij( int, int, double, double**, int );
  double bondorder( int, int, double*, double, double, double**, int );
  double b        ( int, int, double*, double, double, double**, int );

  double gSpline( double, double* );
  double hSpline( double, double* );
  void g_decompose_x( double, size_t*, double* );
  double F_conj( double, double, double, double*, double*, double* );

  void read_file( char * );

  void spline_init();

  void allocate();

  // ----------------------------------------------------------------------
  // S'(t) and S(t) cutoff functions
  // added to header for inlining
  // ----------------------------------------------------------------------

  /* ----------------------------------------------------------------------
     short range cutoff function
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double f_c(double Xij, double Xmin, double Xmax, double *dX) const {
    double cutoff;

    double t = (Xij-Xmin) / (Xmax-Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      *dX = 0.0;
    } else if (t >= 1.0) {
      cutoff = 0.0;
      *dX = 0.0;
    } else {
      double z = t*t*t-1;
      cutoff = exp( gamma_1*t*t*t/z );
      *dX = cutoff * (-3*gamma_1*t*t)/z/z / (Xmax-Xmin);
    }
    return cutoff;
  };

  /* ----------------------------------------------------------------------
     long range cutoff function
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double f_c_LR(double Xij, double Xmin, double Xmax, double *dX) const {
    double cutoff;

    double t = (Xij-Xmin) / (Xmax-Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      //dX = 0.0; this way the derivative is inherited from previous cut off function call
    } else if (t >= 1.0) {
      cutoff = 0.0;
      *dX = 0.0;
    } else {
      cutoff = ( 1.0+cos(MathConst::MY_PI*t) )/2.0;
      *dX = -MathConst::MY_PI*sin(MathConst::MY_PI*t)/2/(Xmax-Xmin);
    }
    return cutoff;
  };

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

UNDOCUMENTED

E: Pair style LCBOP requires atom IDs

UNDOCUMENTED

E: Pair style LCBOP requires newton pair on

UNDOCUMENTED

E: All pair coeffs are not set

UNDOCUMENTED

E: Neighbor list overflow, boost neigh_modify one

UNDOCUMENTED

E: Cannot open LCBOP potential file %s

UNDOCUMENTED

U: Variable name for thermo every does not exist

Self-explanatory.

U: Variable for thermo every is invalid style

Only equal-style variables can be used.

U: Variable name for dump every does not exist

Self-explanatory.

U: Variable for dump every is invalid style

Only equal-style variables can be used.

U: Variable name for restart does not exist

Self-explanatory.

U: Variable for restart is invalid style

Only equal-style variables can be used.

U: Dump every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

U: Restart variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

U: Thermo every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

U: Reuse of dump ID

A dump ID cannot be used twice.

U: Could not find dump group ID

A group ID used in the dump command does not exist.

U: Invalid dump frequency

Dump frequency must be 1 or greater.

U: Invalid dump style

The choice of dump style is unknown.

U: Cound not find dump_modify ID

Self-explanatory.

U: Could not find undump ID

A dump ID used in the undump command does not exist.

U: Thermo_style command before simulation box is defined

The thermo_style command cannot be used before a read_data,
read_restart, or create_box command.

U: New thermo_style command, previous thermo_modify settings will be lost

If a thermo_style command is used after a thermo_modify command, the
settings changed by the thermo_modify command will be reset to their
default values.  This is because the thermo_modify commmand acts on
the currently defined thermo style, and a thermo_style command creates
a new style.

U: Both restart files must have '%' or neither

UNDOCUMENTED

U: Both restart files must use MPI-IO or neither

UNDOCUMENTED

*/
