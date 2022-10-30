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

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifndef LMP_FIX_NH_GPU_H
#define LMP_FIX_NH_GPU_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHGPU : public FixNH {
 public:
  FixNHGPU(class LAMMPS *, int, char **);
  void setup(int vflag) override;
  void reset_dt() override;
  void final_integrate() override;
  double memory_usage() override;

 protected:
  double *_dtfm;
  int _nlocal3, _nlocal_max, _respa_on;

  void remap() override;
  void nve_x() override;
  void nve_v() override;
  void nh_v_press() override;
  void nh_v_temp() override;
};

}    // namespace LAMMPS_NS

#endif
