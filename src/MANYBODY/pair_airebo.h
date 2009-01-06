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

#ifndef PAIR_AIREBO_H
#define PAIR_AIREBO_H

#include "pair.h"

namespace LAMMPS_NS {

class PairAIREBO : public Pair {
 public:
  PairAIREBO(class LAMMPS *);
  ~PairAIREBO();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 private:
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
  void FREBO(int, int);
  void FLJ(int, int);
  void TORSION(int, int);

  double bondorder(int, int, double *, double, double, double **, int);
  double bondorderLJ(int, int, double *, double, double,
		     double *, double, double **, int);

  double Sp(double, double, double, double &);
  double Sp2(double, double, double, double &);

  double gSpline(double, double, int, double *, double *);
  double PijSpline(double, double, int, int, double *);
  double piRCSpline(double, double, double, int, int, double *);
  double TijSpline(double, double, double, double *);

  double kronecker(int, int);

  void add_pages(int);
  void read_file(char *);

  double Sp5th(double, double *, double *);
  double Spbicubic(double, double, double *, double *);
  double Sptricubic(double, double, double, double *, double *);
  void spline_init();

  void allocate();
};

}

#endif
