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

#ifndef MLIAP_DESCRIPTOR_SNAP_H
#define MLIAP_DESCRIPTOR_SNAP_H

#include "pair.h"

namespace LAMMPS_NS {

  //class MLIAPDescriptorSNAP : protected Pointers  {
class MLIAPDescriptorSNAP : public Pair  {
public:
  MLIAPDescriptorSNAP(class LAMMPS *);
  ~MLIAPDescriptorSNAP();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual double memory_usage();

  double rcutfac, quadraticflag; // declared public to workaround gcc 4.9
  int ncoeff;                    //  compiler bug, manifest in KOKKOS package

protected:
  int ncoeffq, ncoeffall;
  class SNA* snaptr;
  virtual void allocate();
  void read_files(char *, char *);
  inline int equal(double* x,double* y);
  inline double dist2(double* x,double* y);

  void compute_beta();
  void compute_bispectrum();

  double rcutmax;               // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  double *radelem;              // element radii
  double *wjelem;               // elements weights
  double **coeffelem;           // element bispectrum coefficients
  double** beta;                // betas for all atoms in list
  double** bispectrum;          // bispectrum components for all atoms in list
  int *map;                     // mapping from atom types to elements
  int twojmax, switchflag, bzeroflag, bnormflag;
  int alloyflag, wselfallflag;
  int chunksize;
  double rfac0, rmin0, wj1, wj2;
  int rcutfacflag, twojmaxflag; // flags for required parameters
  int beta_max;                 // length of beta
};

}

#endif

