/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
NPairStyle(half/respa/bin/newtoff/omp,
           NPairHalfRespaBinNewtoffOmp,
           NP_HALF | NP_RESPA | NP_BIN | NP_NEWTOFF | NP_OMP |
           NP_ORTHO | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_HALF_RESPA_BIN_NEWTOFF_OMP_H
#define LMP_NPAIR_HALF_RESPA_BIN_NEWTOFF_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairHalfRespaBinNewtoffOmp : public NPair {
 public:
  NPairHalfRespaBinNewtoffOmp(class LAMMPS *);
  ~NPairHalfRespaBinNewtoffOmp() {}
  void build(class NeighList *);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

*/
