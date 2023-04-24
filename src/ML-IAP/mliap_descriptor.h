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

#ifndef LMP_MLIAP_DESCRIPTOR_H
#define LMP_MLIAP_DESCRIPTOR_H

#include "pointers.h"

namespace LAMMPS_NS {

class MLIAPDescriptor : virtual protected Pointers {
 public:
  MLIAPDescriptor(LAMMPS *);
  ~MLIAPDescriptor() override;
  virtual void compute_descriptors(class MLIAPData *) = 0;
  virtual void compute_forces(class MLIAPData *) = 0;
  virtual void compute_force_gradients(class MLIAPData *) = 0;
  virtual void compute_descriptor_gradients(class MLIAPData *) = 0;
  virtual void init() = 0;
  virtual double memory_usage();

  int ndescriptors;     // number of descriptors
  int nelements;        // # of unique elements
  char **elements;      // names of unique elements
  double **cutsq;       // nelem x nelem rcutsq values
  double **cutghost;    // cutoff for each ghost pair
  double cutmax;        // maximum cutoff needed
  double *radelem;      // element radii
  double *wjelem;       // elements weights
 protected:
};

}    // namespace LAMMPS_NS

#endif
