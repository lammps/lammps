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
typedef NPairBinOmp<0, 1, 0, 0> NPairFullBinOmp;
NPairStyle(full/bin/omp,
           NPairFullBinOmp,
           NP_FULL | NP_BIN | NP_OMP | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinOmp<1, 0, 0, 0> NPairHalfBinNewtoffOmp;
NPairStyle(half/bin/newtoff/omp,
           NPairHalfBinNewtoffOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinOmp<1, 1, 0, 0> NPairHalfBinNewtonOmp;
NPairStyle(half/bin/newton/omp,
           NPairHalfBinNewtonOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBinOmp<1, 1, 1, 0> NPairHalfBinNewtonTriOmp;
NPairStyle(half/bin/newton/tri/omp,
           NPairHalfBinNewtonTriOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_TRI);

typedef NPairBinOmp<0, 1, 0, 1> NPairFullSizeBinOmp;
NPairStyle(full/size/bin/omp,
           NPairFullSizeBinOmp,
           NP_FULL | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinOmp<1, 0, 0, 1> NPairHalfSizeBinNewtoffOmp;
NPairStyle(half/size/bin/newtoff/omp,
           NPairHalfSizeBinNewtoffOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinOmp<1, 1, 0, 1> NPairHalfSizeBinNewtonOmp;
NPairStyle(half/size/bin/newton/omp,
           NPairHalfSizeBinNewtonOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBinOmp<1, 1, 1, 1> NPairHalfSizeBinNewtonTriOmp;
NPairStyle(half/size/bin/newton/tri/omp,
           NPairHalfSizeBinNewtonTriOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_BIN_OMP_H
#define LMP_NPAIR_BIN_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE>
class NPairBinOmp : public NPair {
 public:
  NPairBinOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
