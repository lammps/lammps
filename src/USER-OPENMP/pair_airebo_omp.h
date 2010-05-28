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

PairStyle(airebo/omp,PairAIREBOOMP)

#else

#ifndef LMP_PAIR_AIREBO_OMP_H
#define LMP_PAIR_AIREBO_OMP_H

#include "pair_omp.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944
#endif

namespace LAMMPS_NS {

class PairAIREBOOMP : public PairOMP {
 public:
  PairAIREBOOMP(class LAMMPS *);
  ~PairAIREBOOMP();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 protected:
  int me;
  int ljflag,torflag;              // 0/1 if LJ,torsion terms included
  int maxlocal;                    // size of numneigh, firstneigh arrays
  int **pages;                     // neighbor list pages
  int maxpage;                     // # of pages currently allocated
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  int npage;                       // current page in page list
  int *map;                        // 0 (C), 1 (H), or -1 (NULL) for each type
  double cutlj;                    // user-specified LJ cutoff
  double cutljrebosq;              // cut for when to compute
                                   // REBO neighs of ghost atoms

  double **cutljsq;                // LJ cutoffs for C,H types
  double **lj1,**lj2,**lj3,**lj4;  // pre-computed LJ coeffs for C,H types
  double cut3rebo;                 // maximum distance for 3rd REBO neigh

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
  double piCCf[5][5][10],piCCdfdx[5][5][10];
  double piCCdfdy[5][5][10],piCCdfdz[5][5][10];
  double piCHf[5][5][10],piCHdfdx[5][5][10];
  double piCHdfdy[5][5][10],piCHdfdz[5][5][10];
  double piHHf[5][5][10],piHHdfdx[5][5][10];
  double piHHdfdy[5][5][10],piHHdfdz[5][5][10];
  double Tf[5][5][10],Tdfdx[5][5][10],Tdfdy[5][5][10],Tdfdz[5][5][10];

  void REBO_neigh();

  template <int EVFLAG, int EFLAG, int VFLAG_ATOM, int NEWTON_PAIR> void FREBO();
  template <int EVFLAG, int EFLAG, int VFLAG_ATOM, int NEWTON_PAIR> void FLJ();
  void TORSION(int, int);

  template <int VFLAG_ATOM>
  double bondorder(int, int, double *, double, double, double **);
  template <int VFLAG_ATOM>
  double bondorderLJ(int, int, double *, double, double,
		     double *, double, double **);

  // ----------------------------------------------------------------------
  // S'(t) and S(t) cutoff functions
  // ----------------------------------------------------------------------

  /* ----------------------------------------------------------------------
     cutoff function Sprime
     return cutoff and dX = derivative
  ------------------------------------------------------------------------- */

  double Sp(double Xij, double Xmin, double Xmax, double &dX) const
    {
      double cutoff;
      
      double t = (Xij-Xmin) / (Xmax-Xmin);
      if (t <= 0.0) {
	cutoff = 1.0;
	dX = 0.0;
      } 
      else if (t >= 1.0) {
	cutoff = 0.0;
	dX = 0.0;
      } 
      else {
	cutoff = 0.5 * (1.0+cos(M_PI*t));
	dX = (-0.5*M_PI*sin(M_PI*t)) / (Xmax-Xmin);
      }
      return cutoff;
    };

/* ----------------------------------------------------------------------
   LJ cutoff function Sp2
   return cutoff and dX = derivative
------------------------------------------------------------------------- */

  double Sp2(double Xij, double Xmin, double Xmax, double &dX) const
    {
      double cutoff;

      double t = (Xij-Xmin) / (Xmax-Xmin);
      if (t <= 0.0) {
	cutoff = 1.0;
	dX = 0.0;
      }
      if (t >= 1.0) {
	cutoff = 0.0;
	dX = 0.0;
      } 
      if (t>0.0 && t<1.0) {
	cutoff = (1.0-(t*t*(3.0-2.0*t)));
	dX = 6.0*(t*t-t) / (Xmax-Xmin);
      }
      return cutoff;
    };

/* ----------------------------------------------------------------------
   fifth order spline evaluation
------------------------------------------------------------------------- */

  double Sp5th(double x, double coeffs[6], double *df) const
    {
      double f;
      int i;
      i = 0;
      f = 0.0;
      *df = 0.0;

      for (i = 0; i<6; i++) {
	f += coeffs[i]*pow(x,((double) i));
	if (i > 0) *df += coeffs[i]*((double) i)*pow(x,((double) i-1.0));
      }

      return f;
    }


  double gSpline(double, double, int, double *, double *);
  double PijSpline(double, double, int, int, double *);
  double piRCSpline(double, double, double, int, int, double *);
  double TijSpline(double, double, double, double *);

/* ----------------------------------------------------------------------
   Kronecker delta function
------------------------------------------------------------------------- */

  double kronecker(const int a, const int b) const
    { return (a == b) ? 1.0 : 0.0; };

  void add_pages(int);
  void read_file(char *);

  double Spbicubic(double, double, double *, double *);
  double Sptricubic(double, double, double, double *, double *);
  void spline_init();

  void allocate();
};

}

#endif
#endif
