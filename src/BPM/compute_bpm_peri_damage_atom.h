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
ComputeStyle(bpm/peri/damage/atom,ComputeBPMPeriDamageAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_BPM_PERI_DAMAGE_ATOM_H
#define LMP_COMPUTE_BPM_PERI_DAMAGE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeBPMPeriDamageAtom : public Compute {
 public:
  ComputeBPMPeriDamageAtom(class LAMMPS *, int, char **);
  ~ComputeBPMPeriDamageAtom() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax;
  int index_vfrac, index_vinter;
  double *damage;
};

}    // namespace LAMMPS_NS

#endif
#endif
