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

#ifdef NPAIR_CLASS
// clang-format off
typedef NPairRespaBin<0, 0> NPairHalfRespaBinNewtoff;
NPairStyle(half/respa/bin/newtoff,
           NPairHalfRespaBinNewtoff,
           NP_HALF | NP_RESPA | NP_BIN | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairRespaBin<1, 0> NPairHalfRespaBinNewton;
NPairStyle(half/respa/bin/newton,
           NPairHalfRespaBinNewton,
           NP_HALF | NP_RESPA | NP_BIN | NP_NEWTON | NP_ORTHO);

typedef NPairRespaBin<1, 1> NPairHalfRespaBinNewtonTri;
NPairStyle(half/respa/bin/newton/tri,
           NPairHalfRespaBinNewtonTri,
           NP_HALF | NP_RESPA | NP_BIN | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_RESPA_BIN_H
#define LMP_NPAIR_RESPA_BIN_H

#include "npair.h"

namespace LAMMPS_NS {

template<int NEWTON, int TRI>
class NPairRespaBin : public NPair {
 public:
  NPairRespaBin(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
