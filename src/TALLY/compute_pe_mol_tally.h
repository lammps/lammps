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
ComputeStyle(pe/mol/tally,ComputePEMolTally);
// clang-format on
#else

#ifndef LMP_COMPUTE_PE_MOL_TALLY_H
#define LMP_COMPUTE_PE_MOL_TALLY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePEMolTally : public Compute {

 public:
  ComputePEMolTally(class LAMMPS *, int, char **);
  ~ComputePEMolTally() override;

  void init() override;
  void compute_vector() override;

  void pair_setup_callback(int, int) override;
  void pair_tally_callback(int, int, int, int, double, double, double, double, double,
                           double) override;

 private:
  bigint did_setup;
  int igroup2, groupbit2;
  double etotal[4];
};

}    // namespace LAMMPS_NS

#endif
#endif
