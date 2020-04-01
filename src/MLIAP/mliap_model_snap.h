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

#ifndef LMP_MLIAP_DESCRIPTOR_SNAP_H
#define LMP_MLIAP_DESCRIPTOR_SNAP_H

#include "pointers.h"

namespace LAMMPS_NS {

class MLIAPModelSNAP : protected Pointers {
public:
  MLIAPModelSNAP(LAMMPS*);
  ~MLIAPModelSNAP();
  virtual void gradient(class NeighList*, double**, double*, double**, int);
  void settings(int, char **);
  virtual void init(int, char **);
  virtual double memory_usage();

  double quadraticflag; // declared public to workaround gcc 4.9
  int ncoeff;           // compiler bug, manifest in KOKKOS package

protected:
  int allocated;
  int ncoeffq, ncoeffall;
  virtual void allocate();
  void read_files(char *, char *);
  inline int equal(double* x,double* y);
  inline double dist2(double* x,double* y);

  void compute_beta(NeighList*);

  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  double **coeffelem;           // element bispectrum coefficients
  double** beta;                // betas for all atoms in list
  double** bispectrum;          // bispectrum components for all atoms in list
  int *map;                     // mapping from atom types to elements
  int twojmax;
  int alloyflag, wselfallflag;
  int twojmaxflag; // flags for required parameters
  int beta_max;                 // length of beta
  int chunksize;
};

}

#endif

