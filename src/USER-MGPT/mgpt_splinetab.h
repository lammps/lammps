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

#ifndef SPLINETAB__
#define SPLINETAB__


/*
  Given a table of ntab data points tab, assumed to be sampled
  on an equidistant grid, compute coefficients of interpolating
  cubic polynimials, one per interval (i.e. ntab-1 polynomials).

  Input point i is located at tab[i*stride].

  Coefficients of output polynomial j are at C[j][0..3];

  The piecewise polynimials form a C^2 function which
  approximates the input function to fourth order.

  The computational cost of this routine is O(ntab).
*/
void makespline(int ntab,int stride,double tab[],double C[][4]);


/*
  Evaluate the spline function with coefficients in C (as returned
  by makespline()) in point x.
  x0 and x1 are the end points of the x points corresponding to
  original input interval tab of makespline(). n is ntab-1.

  The output is the value (y) of the interpolating spline, and the
  first (dy) and second (d2y) dervatives.

  The computational cost of this routine is O(1).
*/
void evalspline(int n,double x0,double x1,double C[][4],
		double x,double *y,double *dy,double *d2y);


/* Evaluate cubic polynomial represented by p in point x.
   The first and second derivatives are also returned.
   This can be used to evaluate one of the sub-polynomials
   in a spline. */
void evalcubic(double p[4],double x,double *y,double *dy,double *d2y);

#endif

