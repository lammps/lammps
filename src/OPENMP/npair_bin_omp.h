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
using NPairFullBinOmp = NPairBinOmp<0, 1, 0, 0, 0>;
NPairStyle(full/bin/omp,
           NPairFullBinOmp,
           NP_FULL | NP_BIN | NP_OMP | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinNewtoffOmp = NPairBinOmp<1, 0, 0, 0, 0>;
NPairStyle(half/bin/newtoff/omp,
           NPairHalfBinNewtoffOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinNewtonOmp = NPairBinOmp<1, 1, 0, 0, 0>;
NPairStyle(half/bin/newton/omp,
           NPairHalfBinNewtonOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfBinNewtonTriOmp = NPairBinOmp<1, 1, 1, 0, 0>;
NPairStyle(half/bin/newton/tri/omp,
           NPairHalfBinNewtonTriOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_TRI);

using NPairFullSizeBinOmp = NPairBinOmp<0, 1, 0, 1, 0>;
NPairStyle(full/size/bin/omp,
           NPairFullSizeBinOmp,
           NP_FULL | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinNewtoffOmp = NPairBinOmp<1, 0, 0, 1, 0>;
NPairStyle(half/size/bin/newtoff/omp,
           NPairHalfSizeBinNewtoffOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinNewtonOmp = NPairBinOmp<1, 1, 0, 1, 0>;
NPairStyle(half/size/bin/newton/omp,
           NPairHalfSizeBinNewtonOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeBinNewtonTriOmp = NPairBinOmp<1, 1, 1, 1, 0>;
NPairStyle(half/size/bin/newton/tri/omp,
           NPairHalfSizeBinNewtonTriOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_MOLONLY | NP_NEWTON | NP_TRI);

using NPairFullBinAtomonlyOmp = NPairBinOmp<0, 1, 0, 0, 1>;
NPairStyle(full/bin/atomonly/omp,
           NPairFullBinAtomonlyOmp,
           NP_FULL | NP_BIN | NP_OMP | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinNewtoffAtomonlyOmp = NPairBinOmp<1, 0, 0, 0, 1>;
NPairStyle(half/bin/atomonly/newtoff/omp,
           NPairHalfBinNewtoffAtomonlyOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinNewtonAtomonlyOmp = NPairBinOmp<1, 1, 0, 0, 1>;
NPairStyle(half/bin/atomonly/newton/omp,
           NPairHalfBinNewtonAtomonlyOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfBinNewtonTriAtomonlyOmp = NPairBinOmp<1, 1, 1, 0, 1>;
NPairStyle(half/bin/newton/tri/atomonly/omp,
           NPairHalfBinNewtonTriAtomonlyOmp,
           NP_HALF | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_TRI);

using NPairFullSizeBinAtomonlyOmp = NPairBinOmp<0, 1, 0, 1, 1>;
NPairStyle(full/size/bin/atomonly/omp,
           NPairFullSizeBinAtomonlyOmp,
           NP_FULL | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinNewtoffAtomonlyOmp = NPairBinOmp<1, 0, 0, 1, 1>;
NPairStyle(half/size/bin/atomonly/newtoff/omp,
           NPairHalfSizeBinNewtoffAtomonlyOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinNewtonAtomonlyOmp = NPairBinOmp<1, 1, 0, 1, 1>;
NPairStyle(half/size/bin/atomonly/newton/omp,
           NPairHalfSizeBinNewtonAtomonlyOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeBinNewtonTriAtomonlyOmp = NPairBinOmp<1, 1, 1, 1, 1>;
NPairStyle(half/size/bin/newton/tri/atomonly/omp,
           NPairHalfSizeBinNewtonTriAtomonlyOmp,
           NP_HALF | NP_SIZE | NP_BIN | NP_OMP | NP_ATOMONLY | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_BIN_OMP_H
#define LMP_NPAIR_BIN_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE, int ATOMONLY>
class NPairBinOmp : public NPair {
 public:
  NPairBinOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
