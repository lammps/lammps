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

/* ----------------------------------------------------------------------
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#include "mgpt_splinetab.h"

static void trisolve(int n,double A[][3],double y[]) {
  /* Backward elimination */
  for(int i = n-1; i>0; i--) {
    double q = A[i-1][2] / A[i][1];
    A[i-1][1] = A[i-1][1] - q*A[i][0];
    y[i-1] = y[i-1] - q*y[i];
  }

  /* Forward substitution */
  y[0] = y[0] / A[0][1];
  for(int i = 1; i<n; i++)
    y[i] = (y[i] - A[i][0]*y[i-1]) / A[i][1];
}

void makespline(int ntab,int stride,double tab[],double C[][4]) {
  int n = 3*(ntab-1);

  double (*A)[3] = new double[n][3];
  double *y = new double[n];

  double h_left,h_right,d;
  int i,j;

  /* Second order second derivative approximation
     at end points. */
  h_left  =
    2.0*tab[stride*0] - 5.0*tab[stride*1] +
    4.0*tab[stride*2] - 1.0*tab[stride*3];
  h_right =
    2.0*tab[stride*(ntab-1)] - 5.0*tab[stride*(ntab-2)] +
    4.0*tab[stride*(ntab-3)] - 1.0*tab[stride*(ntab-4)];

  A[0][0] = 0.0;  A[0][1] = 0.0;  A[0][2] = 2.0;  y[0] = h_left;
  for(i = 1; i<ntab-1; i++) {
    j = 3*(i-1);
    d = tab[stride*i] - tab[stride*(i-1)];
    A[j+1][0] = 1.0;  A[j+1][1] = 1.0;  A[j+1][2] =  1.0;  y[j+1] = d;
    A[j+2][0] = 1.0;  A[j+2][1] = 2.0;  A[j+2][2] = -1.0;  y[j+2] = -d;
    A[j+3][0] = 2.0;  A[j+3][1] = 2.0;  A[j+3][2] = -2.0;  y[j+3] = 2.0*d;
  }

  j = 3*(ntab-2);
  d = tab[stride*(ntab-1)] - tab[stride*(ntab-2)];
  A[j+1][0] = 1.0;  A[j+1][1] = 1.0;  A[j+1][2] =  1.0;  y[j+1] = d;
  A[j+2][0] = 2.0;  A[j+2][1] = 6.0;  A[j+2][2] =  0.0;  y[j+2] = h_right;

  trisolve(n,A,y);

  for(i = 0; i<ntab-1; i++) {
    C[i][0] = tab[stride*i];
    C[i][1] = y[3*i+0];
    C[i][2] = y[3*i+1];
    C[i][3] = y[3*i+2];
  }

  delete[] y;
  delete[] A;
}

void evalcubic(double p[4],double x,double *y,double *dy,double *d2y) {
  double t1,t2,t3;

  t1 = p[2] + x*p[3];
  t2 = p[1] + x*t1;

  t3 = t1 + x*p[3];

  *y = p[0] + x*t2;
  *dy = (t2 + x*t3);
  *d2y = 2.0*(t3 + x*p[3]);
}

void evalspline(int n,double x0,double x1,double C[][4],
		  double x,double *y,double *dy,double *d2y) {
  double xhat,t1,t2,t3;
  double *p;
  int idx;
  double dxinv = n/(x1-x0);
  xhat = (x-x0)/(x1-x0) * n;

  idx = (int) xhat;
  if(idx < 0) idx = 0;
  if(idx > n-1) idx = n-1;
  xhat = xhat - idx;
  p = C[idx];

  if(0) {
    *y = p[0] + xhat*(p[1] + xhat*(p[2] + xhat*p[3]));

    *dy = p[1] + xhat*(2*p[2] + xhat*3*p[3]);
    *d2y = 2*p[2] + xhat*6*p[3];

    *dy *= dxinv;
    *d2y *= dxinv*dxinv;
  } else {
    t1 = p[2] + xhat*p[3];
    t2 = p[1] + xhat*t1;

    t3 = t1 + xhat*p[3];

    *y = p[0] + xhat*t2;
    *dy = (t2 + xhat*t3)*dxinv;
    *d2y = 2.0*(t3 + xhat*p[3])*(dxinv*dxinv);
  }
}
