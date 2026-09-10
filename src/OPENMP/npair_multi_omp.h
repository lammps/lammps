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
using NPairFullMultiOmp = NPairMultiOmp<0, 1, 0, 0, 0>;
NPairStyle(full/multi/omp,
           NPairFullMultiOmp,
           NP_FULL | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiNewtoffOmp = NPairMultiOmp<1, 0, 0, 0, 0>;
NPairStyle(half/multi/newtoff/omp,
           NPairHalfMultiNewtoffOmp,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiNewtonOmp = NPairMultiOmp<1, 1, 0, 0, 0>;
NPairStyle(half/multi/newton/omp,
           NPairHalfMultiNewtonOmp,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTON | NP_ORTHO);

using NPairHalfMultiNewtonTriOmp = NPairMultiOmp<1, 1, 1, 0, 0>;
NPairStyle(half/multi/newton/tri/omp,
           NPairHalfMultiNewtonTriOmp,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTON | NP_TRI);

using NPairFullSizeMultiOmp = NPairMultiOmp<0, 1, 0, 1, 0>;
NPairStyle(full/size/multi/omp,
           NPairFullSizeMultiOmp,
           NP_FULL | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiNewtoffOmp = NPairMultiOmp<1, 0, 0, 1, 0>;
NPairStyle(half/size/multi/newtoff/omp,
           NPairHalfSizeMultiNewtoffOmp,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiNewtonOmp = NPairMultiOmp<1, 1, 0, 1, 0>;
NPairStyle(half/size/multi/newton/omp,
           NPairHalfSizeMultiNewtonOmp,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeMultiNewtonTriOmp = NPairMultiOmp<1, 1, 1, 1, 0>;
NPairStyle(half/size/multi/newton/tri/omp,
           NPairHalfSizeMultiNewtonTriOmp,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_OMP | NP_NEWTON | NP_TRI);

using NPairFullMultiAtomonlyOmp = NPairMultiOmp<0, 1, 0, 0, 1>;
NPairStyle(full/multi/atomonly/omp,
           NPairFullMultiAtomonlyOmp,
           NP_FULL | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiAtomonlyNewtoffOmp = NPairMultiOmp<1, 0, 0, 0, 1>;
NPairStyle(half/multi/atomonly/newtoff/omp,
           NPairHalfMultiAtomonlyNewtoffOmp,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiAtomonlyNewtonOmp = NPairMultiOmp<1, 1, 0, 0, 1>;
NPairStyle(half/multi/atomonly/newton/omp,
           NPairHalfMultiAtomonlyNewtonOmp,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTON | NP_ORTHO);

using NPairHalfMultiAtomonlyNewtonTriOmp = NPairMultiOmp<1, 1, 1, 0, 1>;
NPairStyle(half/multi/atomonly/newton/tri/omp,
           NPairHalfMultiAtomonlyNewtonTriOmp,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTON | NP_TRI);

using NPairFullSizeMultiAtomonlyOmp = NPairMultiOmp<0, 1, 0, 1, 1>;
NPairStyle(full/size/multi/atomonly/omp,
           NPairFullSizeMultiAtomonlyOmp,
           NP_FULL | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiAtomonlyNewtoffOmp = NPairMultiOmp<1, 0, 0, 1, 1>;
NPairStyle(half/size/multi/atomonly/newtoff/omp,
           NPairHalfSizeMultiAtomonlyNewtoffOmp,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiAtomonlyNewtonOmp = NPairMultiOmp<1, 1, 0, 1, 1>;
NPairStyle(half/size/multi/atomonly/newton/omp,
           NPairHalfSizeMultiAtomonlyNewtonOmp,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeMultiAtomonlyNewtonTriOmp = NPairMultiOmp<1, 1, 1, 1, 1>;
NPairStyle(half/size/multi/atomonly/newton/tri/omp,
           NPairHalfSizeMultiAtomonlyNewtonTriOmp,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_OMP | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_MULTI_OMP_H
#define LMP_NPAIR_MULTI_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE, int ATOMONLY>
class NPairMultiOmp : public NPair {
 public:
  NPairMultiOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
