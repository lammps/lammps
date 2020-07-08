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

#ifdef COMPUTE_CLASS

ComputeStyle(mliap,ComputeMLIAP)

#else

#ifndef LMP_COMPUTE_MLIAP_H
#define LMP_COMPUTE_MLIAP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMLIAP : public Compute {
 public:
  ComputeMLIAP(class LAMMPS *, int, char **);
  ~ComputeMLIAP();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  void compute_array_alt();
  double memory_usage();

 private:
  int natoms, nmax, size_gradforce, lastcol;
  int yoffset, zoffset;
  int ndims_force, ndims_virial;
  class NeighList *list;
  double **mliap, **mliapall;
  double **gradforce;
  int *map;  // map types to [0,nelements)
  int nelements;
  int gradgradflag;           // 1 for graddesc, 0 for gamma
  double** descriptors;        // descriptors for all atoms in list
  int ndescriptors;            // number of descriptors 
  int nparams;                 // number of model parameters per element
  int gamma_nnz;               // number of non-zero entries in gamma
  double** gamma;              // gamma element
  int** gamma_row_index;       // row (parameter) index 
  int** gamma_col_index;       // column (descriptor) index 
  double* egradient;           // energy gradient w.r.t. parameters

  // data structures for descriptor neighbor list
  // this is for neighbors strictly inside deescriptor cutoff

  int numlistdesc_max;         // number of atoms allocated in neighlist
  int *numneighdesc;           // numbers of neighbors for descriptors
  int *iatomdesc;              // list of descriptor atoms
  int *ielemdesc;              // list of descriptor elements
  int nneighdesc_max;          // number of ij neighbors allocated in graddesc
  int *jatomdesc;              // list of descriptor neighbor atoms
  int *jelemdesc;              // list of descriptor neighbor elements
  int *neighdesc;              // list of descriptor neighbors
  double ***graddesc;          // gradient of descriptors w.r.t. rij

  class MLIAPModel* model;
  class MLIAPDescriptor* descriptor;

  Compute *c_pe;
  Compute *c_virial;

  void dbdotr_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute snap requires a pair style be defined

Self-explanatory.

E: Compute snap cutoff is longer than pairwise cutoff

UNDOCUMENTED

W: More than one compute snad/atom

Self-explanatory.

*/
