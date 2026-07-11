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

#ifdef INTEGRATE_CLASS
// clang-format off
IntegrateStyle(respa/omp,RespaOMP);
// clang-format on
#else

#ifndef LMP_RESPA_OMP_H
#define LMP_RESPA_OMP_H

#include "respa.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class RespaOMP : public Respa, public ThrOMP {
 public:
  RespaOMP(class LAMMPS *, int, char **);

  void init() override;
  void setup(int) override;
  void setup_minimal(int) override;

 protected:
  void recurse(int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
