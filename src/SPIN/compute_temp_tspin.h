/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/tspin,ComputeTempTSpin);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_TSPIN_H
#define LMP_COMPUTE_TEMP_TSPIN_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempTSpin : public Compute {
 public:
  ComputeTempTSpin(class LAMMPS *, int, char **);

  void init() override;
  double compute_scalar() override;

 private:
  int index_vs, index_sm;
  double tfactor;
};

}    // namespace LAMMPS_NS

#endif
#endif
