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
typedef NPairBinAtomonlyOmp<0, 1, 0, 0> NPairFullBinAtomonlyOmp;
NPairStyle(full/bin/atomonly/omp,
           NPairFullBinAtomonlyOmp,
           NP_FULL | NP_BIN | NP_OMP | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinAtomonlyOmp<1, 0, 0, 0> NPairHalfBinNewtoffAtomonlyOmp;
NPairStyle(half/bin/newtoff/atomonly/omp,
           NPairHalfBinNewtoffAtomonlyOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinAtomonlyOmp<1, 1, 0, 0> NPairHalfBinNewtonAtomonlyOmp;
NPairStyle(half/bin/newton/atomonly/omp,
           NPairHalfBinNewtonAtomonlyOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBinAtomonlyOmp<1, 1, 1, 0> NPairHalfBinNewtonTriAtomonlyOmp;
NPairStyle(half/bin/newton/tri/atomonly/omp,
           NPairHalfBinNewtonTriAtomonlyOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_TRI);

typedef NPairBinAtomonlyOmp<0, 1, 0, 1> NPairFullSizeBinAtomonlyOmp;
NPairStyle(full/size/bin/atomonly/omp,
           NPairFullSizeBinAtomonlyOmp,
           NP_FULL | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinAtomonlyOmp<1, 0, 0, 1> NPairHalfSizeBinNewtoffAtomonlyOmp;
NPairStyle(half/size/bin/newtoff/atomonly/omp,
           NPairHalfSizeBinNewtoffAtomonlyOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinAtomonlyOmp<1, 1, 0, 1> NPairHalfSizeBinNewtonAtomonlyOmp;
NPairStyle(half/size/bin/newton/atomonly/omp,
           NPairHalfSizeBinNewtonAtomonlyOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBinAtomonlyOmp<1, 1, 1, 1> NPairHalfSizeBinNewtonTriAtomonlyOmp;
NPairStyle(half/size/bin/newton/tri/atomonly/omp,
           NPairHalfSizeBinNewtonTriAtomonlyOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_TRI);
// clang-format on
#else


#ifndef LMP_NPAIR_BIN_ATOMONLY_OMP_H
#define LMP_NPAIR_BIN_ATOMONLY_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE>
class NPairBinAtomonlyOmp : public NPair {
 public:
  NPairBinAtomonlyOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
