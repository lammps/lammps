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

#ifndef LMP_MLIAP_H
#define LMP_MLIAP_H

#include "pointers.h"

namespace LAMMPS_NS {

class MLIAP : protected Pointers {

 public:
  MLIAP(class LAMMPS*, int, int, int, int, int*, class MLIAPModel*, class MLIAPDescriptor*);
  ~MLIAP();

  void init();
  void generate_neigharrays(class NeighList *);
  void grow_neigharrays();
  double memory_usage();

  // private:
  int size_array_rows, size_array_cols;
  int natoms, size_gradforce;
  int yoffset, zoffset;
  int ndims_force, ndims_virial;
  class NeighList *list;
  double **gradforce;

  int *map;  // map types to [0,nelements)
  double** beta;                // betas for all atoms in list
  double** descriptors;        // descriptors for all atoms in list
  int ndescriptors;            // number of descriptors 
  int nparams;                 // number of model parameters per element
  int nelements;               // number of elements
  int gradgradflag;            // 1 for graddesc, 0 for gamma
        
  // data structures for grad-grad list (gamma)

  int gamma_nnz;               // number of non-zero entries in gamma
  double** gamma;              // gamma element
  int** gamma_row_index;       // row (parameter) index 
  int** gamma_col_index;       // column (descriptor) index 
  double* egradient;           // energy gradient w.r.t. parameters

  // data structures for mliap neighbor list
  // only neighbors strictly inside descriptor cutoff

  int natomdesc;                // current number of atoms
  int natomdesc_max;            // allocated size of descriptor array
  int natomneigh_max;           // allocated size of atom neighbor arrays
  int *numneighmliap;           // neighbors count for each atom
  int *iatommliap;              // index of each atom
  int *ielemmliap;              // element of each atom
  int nneigh_max;               // number of ij neighbors allocated
  int *jatommliap;              // index of each neighbor
  int *jelemmliap;              // element of each neighbor
  double ***graddesc;           // descriptor gradient w.r.t. each neighbor

 private:
  class MLIAPModel* model;
  class MLIAPDescriptor* descriptor;

};

}

#endif
