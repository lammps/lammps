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
ComputeStyle(rmsd,ComputeRmsd);
// clang-format on
#else

#ifndef LMP_COMPUTE_RMSD_H
#define LMP_COMPUTE_RMSD_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRmsd : public Compute {
 public:
  ComputeRmsd(class LAMMPS *, int, char **);
  ~ComputeRmsd() override;
  void init() override;
  double compute_scalar() override;

 protected:

  tagint *group_taglist;
  double **x_group, **x_group_shifted, **ref_positions_shifted;

  static int idcompare(const tagint, const tagint, void *);

  // -------- RMSD --------

  tagint group_count;
  double **ref_positions;

  double rmsd(double *);

  // -------- PRIVATE IMPLEMENTATION METHODS --------

  void read_xyz(char *);

};
}    // namespace LAMMPS_NS


#endif
#endif
