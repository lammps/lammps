/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MLIAPDATA_H
#define LMP_MLIAPDATA_H

#include "pointers.h"

namespace LAMMPS_NS {

class MLIAPData : protected Pointers {

 public:
  MLIAPData(class LAMMPS *, int, int *, class MLIAPModel *, class MLIAPDescriptor *,
            class PairMLIAP * = nullptr);
  ~MLIAPData() override;

  void init();
  virtual void generate_neighdata(class NeighList *, int = 0, int = 0);
  virtual void grow_neigharrays();
  double memory_usage();

  int size_array_rows, size_array_cols;
  int natoms;
  int size_gradforce;
  int yoffset, zoffset;
  int ndims_force, ndims_virial;
  double **f;
  double **gradforce;
  double **betas;          // betas for all atoms in list
  double **descriptors;    // descriptors for all atoms in list
  double *eatoms;          // energies for all atoms in list
  double energy;           // energy
  int ndescriptors;        // number of descriptors
  int nparams;             // number of model parameters per element
  int nelements;           // number of elements
  int gradgradflag;        // 1 for graddesc, 0 for gamma, -1 for pair style

  // data structures for grad-grad list (gamma)

  int natomgamma_max;       // allocated size of gamma
  int gamma_nnz;            // number of non-zero entries in gamma
  double **gamma;           // gamma element
  int **gamma_row_index;    // row (parameter) index
  int **gamma_col_index;    // column (descriptor) index
  double *egradient;        // energy gradient w.r.t. parameters

  // data structures for mliap neighbor list
  // only neighbors strictly inside descriptor cutoff

  int ntotal;        // total number of owned and ghost atoms on this proc
  int nlistatoms;    // current number of atoms in local atom lists
  int nlocal;
  int nlistatoms_max;            // allocated size of descriptor array
  int natomneigh;                // current number of atoms and ghosts in atom neighbor arrays
  int natomneigh_max;            // allocated size of atom neighbor arrays
  int *numneighs;                // neighbors count for each atom
  int *iatoms;                   // index of each atom
  int *ielems;                   // element of each atom
  int *itypes;                   // LAMMPS type of each atom for external evaluators
  int nneigh_max;                // number of ij neighbors allocated
  int npairs;                    // number of ij neighbor pairs
  int *pair_i;                   // index of each i atom for each ij pair
  int *jatoms;                   // index of each neighbor
  int *jelems;                   // element of each neighbor
  int *elems;                    // element of each atom in or not in the neighborlist
  int **lmp_firstneigh;          // copy of list->firstneigh for external evaluators
  double **rij;                  // distance vector of each neighbor
  double ***graddesc;            // descriptor gradient w.r.t. each neighbor
  int eflag;                     // indicates if energy is needed
  int vflag;                     // indicates if virial is needed
  class PairMLIAP *pairmliap;    // access to pair tally functions

 protected:
  class MLIAPModel *model;
  class MLIAPDescriptor *descriptor;

  int nmax;
  class NeighList *list;    // LAMMPS neighbor list
  int *map;                 // map LAMMPS types to [0,nelements)
};

}    // namespace LAMMPS_NS

#endif
