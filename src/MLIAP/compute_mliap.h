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
  double memory_usage();

 private:
  int natoms, nmax, size_peratom, lastcol;
  int nperdim, yoffset, zoffset;
  int ndims_peratom, ndims_force, ndims_virial;
  double **cutsq;
  class NeighList *list;
  double **mliap, **mliapall;
  double **mliap_peratom;
  int *map;  // map types to [0,nelements)
  int nelements;

  double*** gamma;             // gammas for all atoms in list
  double** descriptors;        // descriptors for all atoms in list
  int ndescriptors;            // number of descriptors 
  int gamma_max;               // number of atoms allocated for beta, descriptors
  int nparams;                 // number of model paramters per element

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
