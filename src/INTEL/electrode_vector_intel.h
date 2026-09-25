/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Shern Tee (UQ)
 ------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_VECTOR_INTEL_H
#define LMP_ELECTRODE_VECTOR_INTEL_H

#include "electrode_vector.h"
#include "fix_intel.h"

namespace LAMMPS_NS {

class ElectrodeVectorIntel : public ElectrodeVector {
 public:
  ElectrodeVectorIntel(class LAMMPS *lmp, int narg, char **arg, int sensor_group,
     int source_group, double eta, bool invert_source) :
  ElectrodeVector(lmp, narg, arg, sensor_group, source_group, eta, invert_source),
  fix(nullptr)
  {
  }

 private:
  void get_fix_intel() override;
  void pair_contribution(double *) override;
  FixIntel *fix;
  int _cop, _lrt, _ccache_stride;
  template <class flt_t, class acc_t> void pair_contribution(IntelBuffers<flt_t, acc_t> *buffers, double *vec);
};

}

#endif
