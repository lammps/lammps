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

#ifndef LMP_MLIAP_DESCRIPTOR_SNAP_H
#define LMP_MLIAP_DESCRIPTOR_SNAP_H

#include "mliap_descriptor.h"

namespace LAMMPS_NS {

class MLIAPDescriptorSNAP : public MLIAPDescriptor {
 public:
  MLIAPDescriptorSNAP(LAMMPS *, char *);
  ~MLIAPDescriptorSNAP() override;
  void compute_descriptors(class MLIAPData *) override;
  void compute_forces(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  void compute_descriptor_gradients(class MLIAPData *) override;
  void init() override;
  double memory_usage() override;

  double rcutfac;

 protected:
  class SNA *snaptr;
  void read_paramfile(char *);
  inline int equal(double *x, double *y);
  inline double dist2(double *x, double *y);

  int twojmax, switchflag, bzeroflag;
  int chemflag, bnormflag, wselfallflag;
  int switchinnerflag;
  double rfac0, rmin0;

  double *sinnerelem;
  double *dinnerelem;
};

}    // namespace LAMMPS_NS

#endif
