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

#ifdef PAIR_CLASS

PairStyle(airebo,PairAIREBO)

#else

#ifndef LMP_PAIR_AIREBO_H
#define LMP_PAIR_AIREBO_H

#include "pair.h"
#include "my_page.h"
#include "math.h"
#include "math_const.h"

namespace LAMMPS_NS {

class PairAIREBO : public Pair {
 public:
  PairAIREBO(class LAMMPS *);
  virtual ~PairAIREBO();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 protected:
  int *map;                        // 0 (C), 1 (H), or -1 (NULL) for each type

  int me;
  int ljflag,torflag;              // 0/1 if LJ,torsion terms included

  double cutlj;                    // user-specified LJ cutoff
  double cutljrebosq;              // cut for when to compute
                                   // REBO neighs of ghost atoms

  double **cutljsq;                // LJ cutoffs for C,H types
  double **lj1,**lj2,**lj3,**lj4;  // pre-computed LJ coeffs for C,H types
  double cut3rebo;                 // maximum distance for 3rd REBO neigh

  int maxlocal;                    // size of numneigh, firstneigh arrays
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  MyPage<int> *ipage;              // neighbor list pages
  int *REBO_numneigh;              // # of pair neighbors for each atom
  int **REBO_firstneigh;           // ptr to 1st neighbor of each atom

  double *closestdistsq;           // closest owned atom dist to each ghost
  double *nC,*nH;                  // sum of weighting fns with REBO neighs

  double smin,Nmin,Nmax,NCmin,NCmax,thmin,thmax;
  double rcmin[2][2],rcmax[2][2],rcmaxsq[2][2],rcmaxp[2][2];
  double Q[2][2],alpha[2][2],A[2][2],rho[2][2],BIJc[2][2][3],Beta[2][2][3];
  double rcLJmin[2][2],rcLJmax[2][2],rcLJmaxsq[2][2],bLJmin[2][2],bLJmax[2][2];
  double epsilon[2][2],sigma[2][2],epsilonT[2][2];

  // spline coefficients

  double gCdom[5],gC1[4][6],gC2[4][6],gHdom[4],gH[3][6];
  double pCCdom[2][2],pCHdom[2][2],pCC[4][4][16],pCH[4][4][16];
  double piCCdom[3][2],piCHdom[3][2],piHHdom[3][2];
  double piCC[4][4][9][64],piCH[4][4][9][64],piHH[4][4][9][64];
  double Tijdom[3][2],Tijc[4][4][9][64];

  // spline knot values

  double PCCf[5][5],PCCdfdx[5][5],PCCdfdy[5][5],PCHf[5][5];
  double PCHdfdx[5][5],PCHdfdy[5][5];
  double piCCf[5][5][11],piCCdfdx[5][5][11];
  double piCCdfdy[5][5][11],piCCdfdz[5][5][11];
  double piCHf[5][5][11],piCHdfdx[5][5][11];
  double piCHdfdy[5][5][11],piCHdfdz[5][5][11];
  double piHHf[5][5][11],piHHdfdx[5][5][11];
  double piHHdfdy[5][5][11],piHHdfdz[5][5][11];
  double Tf[5][5][10],Tdfdx[5][5][10],Tdfdy[5][5][10],Tdfdz[5][5][10];

  void REBO_neigh();
  void FREBO(int, int);
  void FLJ(int, int);
  void TORSION(int, int);

  double bondorder(int, int, double *, double, double, double **, int);
  double bondorderLJ(int, int, double *, double, double,
                     double *, double, double **, int);

  double gSpline(double, double, int, double *, double *);
  double PijSpline(double, double, int, int, double *);
  double piRCSpline(double, double, double, int, int, double *);
  double TijSpline(double, double, double, double *);

  void read_file(char *);

  double Sp5th(double, double *, double *);
  double Spbicubic(double, double, double *, double *);
  double Sptricubic(double, double, double, double *, double *);
  void spline_init();

  void allocate();

  // ----------------------------------------------------------------------
  // S'(t) and S(t) cutoff functions
  // added to header for inlining
  // ----------------------------------------------------------------------

  /* ----------------------------------------------------------------------
     cutoff function Sprime
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double Sp(double Xij, double Xmin, double Xmax, double &dX) const {
    double cutoff;

    double t = (Xij-Xmin) / (Xmax-Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      dX = 0.0;
    } else if (t >= 1.0) {
      cutoff = 0.0;
      dX = 0.0;
    } else {
      cutoff = 0.5 * (1.0+cos(t*MathConst::MY_PI));
      dX = (-0.5*MathConst::MY_PI*sin(t*MathConst::MY_PI)) / (Xmax-Xmin);
    }
    return cutoff;
  };

  /* ----------------------------------------------------------------------
     LJ cutoff function Sp2
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double Sp2(double Xij, double Xmin, double Xmax, double &dX) const {
    double cutoff;

    double t = (Xij-Xmin) / (Xmax-Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      dX = 0.0;
    } else if (t >= 1.0) {
      cutoff = 0.0;
      dX = 0.0;
    } else {
      cutoff = (1.0-(t*t*(3.0-2.0*t)));
      dX = 6.0*(t*t-t) / (Xmax-Xmin);
    }
    return cutoff;
  };

  /* kronecker delta function returning a double */

  inline double kronecker(const int a, const int b) const {
    return (a == b) ? 1.0 : 0.0;
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

Self-explanatory.  Check the input script or data file.

E: Pair style AIREBO requires atom IDs

This is a requirement to use the AIREBO potential.

E: Pair style AIREBO requires newton pair on

See the newton command.  This is a restriction to use the AIREBO
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Neighbor list overflow, boost neigh_modify one

There are too many neighbors of a single atom.  Use the neigh_modify
command to increase the max number of neighbors allowed for one atom.
You may also want to boost the page size.

E: Cannot open AIREBO potential file %s

The specified AIREBO potential file cannot be opened.  Check that the
path and name are correct.

*/
