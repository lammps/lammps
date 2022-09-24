/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_UNIFIED_H
#define LMP_MLIAP_UNIFIED_H

#include "mliap_data.h"
#include "mliap_descriptor.h"
#include "mliap_model.h"

#include <Python.h>

namespace LAMMPS_NS {

class MLIAPDummyDescriptor : public MLIAPDescriptor {
 public:
  MLIAPDummyDescriptor(LAMMPS *);
  ~MLIAPDummyDescriptor();
  virtual void compute_descriptors(class MLIAPData *);
  virtual void compute_forces(class MLIAPData *);
  virtual void compute_force_gradients(class MLIAPData *);
  virtual void compute_descriptor_gradients(class MLIAPData *);
  virtual void init();
  void set_elements(char **, int);

  PyObject *unified_interface;    // MLIAPUnifiedInterface
  double rcutfac;
};

class MLIAPDummyModel : public MLIAPModel {
 public:
  MLIAPDummyModel(LAMMPS *, char * = NULL);
  ~MLIAPDummyModel();
  virtual int get_nparams();
  virtual int get_gamma_nnz(class MLIAPData *);
  virtual void compute_gradients(class MLIAPData *);
  virtual void compute_gradgrads(class MLIAPData *);
  virtual void compute_force_gradients(class MLIAPData *);
  virtual double memory_usage();

  PyObject *unified_interface;

 protected:
  virtual void read_coeffs(char *);
};

struct MLIAPBuildUnified_t {
  MLIAPData *data;
  MLIAPDummyDescriptor *descriptor;
  MLIAPDummyModel *model;
};

MLIAPBuildUnified_t build_unified(char *, MLIAPData *, LAMMPS *, char * = NULL);
void update_pair_energy(MLIAPData *, double *);
void update_pair_forces(MLIAPData *, double *);

}    // namespace LAMMPS_NS

#endif
