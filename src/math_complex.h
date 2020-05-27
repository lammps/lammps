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

/* ----------------------------------------------------------------------
   Contributing author: Pieter J. in 't Veld (SNL)
------------------------------------------------------------------------- */

#ifndef LMP_MATH_COMPLEX_H
#define LMP_MATH_COMPLEX_H

#define COMPLEX_NULL        {0, 0}

namespace LAMMPS_NS {

typedef struct complex {
  double re, im; } complex;

}

#define C_MULT(d, x, y) { \
  d.re = x.re*y.re-x.im*y.im; \
  d.im = x.re*y.im+x.im*y.re; }

#define C_RMULT(d, x, y) { \
  complex t = x; \
  d.re = t.re*y.re-t.im*y.im; \
  d.im = t.re*y.im+t.im*y.re; }

#define C_CRMULT(d, x, y) { \
  complex t = x; \
  d.re = t.re*y.re-t.im*y.im; \
  d.im = -t.re*y.im-t.im*y.re; }

#define C_SMULT(d, x, y) { \
  d.re = x.re*y; \
  d.im = x.im*y; }

#define C_ADD(d, x, y) { \
  d.re = x.re+y.re; \
  d.im = x.im+y.im; }

#define C_SUBTR(d, x, y) { \
  d.re = x.re-y.re; \
  d.im = x.im-y.im; }

#define C_CONJ(d, x) { \
  d.re = x.re; \
  d.im = -x.im; }

#define C_SET(d, x, y) { \
  d.re = x; \
  d.im = y; }

#define C_ANGLE(d, angle) { \
  double a = angle; \
  d.re = cos(a); \
  d.im = sin(a); }

#define C_COPY(d, x) { \
  memcpy(&d, &x, sizeof(complex)); }

#endif
