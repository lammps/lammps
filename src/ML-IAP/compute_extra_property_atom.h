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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(extraProperty/atom,ComputeExtraPropertyAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_EXTRA_PROPERTY_ATOM_H
#define LMP_COMPUTE_EXTRA_PROPERTY_ATOM_H

#include "compute.h"
#include "mliap_data.h"
#include "mliap_descriptor.h"
#include "pair_mliap.h"

namespace LAMMPS_NS {

class ComputeExtraPropertyAtom : public Compute {
 public:
  ComputeExtraPropertyAtom(class LAMMPS *, int, char **);
  ~ComputeExtraPropertyAtom() override;
  void init() override;
  void compute_peratom() override;

 private:
  std::string extra_property_name;
  int hybridIndex;
  MLIAPData *data;
};

}    // namespace LAMMPS_NS

#endif
#endif
