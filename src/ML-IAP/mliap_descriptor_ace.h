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

#ifndef LMP_MLIAP_DESCRIPTOR_ACE_H
#define LMP_MLIAP_DESCRIPTOR_ACE_H

#include "mliap_descriptor.h"

namespace LAMMPS_NS {

class MLIAPDescriptorACE : public MLIAPDescriptor {
 public:
  MLIAPDescriptorACE(LAMMPS *, char *);
  ~MLIAPDescriptorACE() override;
  void compute_descriptors(class MLIAPData *) override;
  void compute_forces(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  void compute_descriptor_gradients(class MLIAPData *) override;
  void init() override;
  double memory_usage() override;

  double rcutfac;
  int allocated = 0;
  int max_num = 0;
  char *ctilde_file;

 protected:
  virtual void allocate();
  int natoms, nmax, size_peratom, lastcol;
  int ncoeff, nvalues, nperdim, yoffset, zoffset;
  int ndims_peratom, ndims_force, ndims_virial;
  int n_r1, n_rp;
  int chemflag;
  int bikflag, bik_rows, dgradflag, dgrad_rows;
  double cutmax;
  struct ACE_ML_impl *acemlimpl;
};

}    // namespace LAMMPS_NS

#endif
