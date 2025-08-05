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
ComputeStyle(mliap/property/atom,ComputeMLIAPPropertyAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_MLIAP_PROPERTY_ATOM_H
#define LMP_COMPUTE_MLIAP_PROPERTY_ATOM_H

#include "compute.h"
#include "mliap_data.h"
#include "mliap_descriptor.h"
#include "pair_mliap.h"

namespace LAMMPS_NS {

class ComputeMLIAPPropertyAtom : public Compute {
 public:
  ComputeMLIAPPropertyAtom(class LAMMPS *, int, char **);
  ~ComputeMLIAPPropertyAtom() override;
  void init() override;
  void compute_peratom() override;

 private:
  std::string property_name;
  int hybrid_index;
  MLIAPData *data;
};

}    // namespace LAMMPS_NS

#endif
#endif
