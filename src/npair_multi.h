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
typedef NPairMulti<0, 1, 0, 0, 0> NPairFullMulti;
NPairStyle(full/multi,
           NPairFullMulti,
           NP_FULL | NP_MULTI | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 0, 0, 0, 0> NPairHalfMultiNewtoff;
NPairStyle(half/multi/newtoff,
           NPairHalfMultiNewtoff,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 1, 0, 0, 0> NPairHalfMultiNewton;
NPairStyle(half/multi/newton,
           NPairHalfMultiNewton,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

typedef NPairMulti<1, 1, 1, 0, 0> NPairHalfMultiNewtonTri;
NPairStyle(half/multi/newton/tri,
           NPairHalfMultiNewtonTri,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_TRI);

typedef NPairMulti<0, 1, 0, 1, 0> NPairFullSizeMulti;
NPairStyle(full/size/multi,
           NPairFullSizeMulti,
           NP_FULL | NP_SIZE | NP_MULTI | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 0, 0, 1, 0> NPairHalfSizeMultiNewtoff;
NPairStyle(half/size/multi/newtoff,
           NPairHalfSizeMultiNewtoff,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 1, 0, 1, 0> NPairHalfSizeMultiNewton;
NPairStyle(half/size/multi/newton,
           NPairHalfSizeMultiNewton,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

typedef NPairMulti<1, 1, 1, 1, 0> NPairHalfSizeMultiNewtonTri;
NPairStyle(half/size/multi/newton/tri,
           NPairHalfSizeMultiNewtonTri,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_TRI);

typedef NPairMulti<0, 1, 0, 0, 1> NPairFullMultiAtomonly;
NPairStyle(full/multi/atomonly,
           NPairFullMultiAtomonly,
           NP_FULL | NP_MULTI | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 0, 0, 0, 1> NPairHalfMultiAtomonlyNewtoff;
NPairStyle(half/multi/atomonly/newtoff,
           NPairHalfMultiAtomonlyNewtoff,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 1, 0, 0, 1> NPairHalfMultiAtomonlyNewton;
NPairStyle(half/multi/atomonly/newton,
           NPairHalfMultiAtomonlyNewton,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

typedef NPairMulti<1, 1, 1, 0, 1> NPairHalfMultiAtomonlyNewtonTri;
NPairStyle(half/multi/atomonly/newton/tri,
           NPairHalfMultiAtomonlyNewtonTri,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_NEWTON | NP_TRI);

typedef NPairMulti<0, 1, 0, 1, 1> NPairFullSizeMultiAtomonly;
NPairStyle(full/size/multi/atomonly,
           NPairFullSizeMultiAtomonly,
           NP_FULL | NP_SIZE | NP_MULTI | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 0, 0, 1, 1> NPairHalfSizeMultiAtomonlyNewtoff;
NPairStyle(half/size/multi/atomonly/newtoff,
           NPairHalfSizeMultiAtomonlyNewtoff,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMulti<1, 1, 0, 1, 1> NPairHalfSizeMultiAtomonlyNewton;
NPairStyle(half/size/multi/atomonly/newton,
           NPairHalfSizeMultiAtomonlyNewton,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

typedef NPairMulti<1, 1, 1, 1, 1> NPairHalfSizeMultiAtomonlyNewtonTri;
NPairStyle(half/size/multi/atomonly/newton/tri,
           NPairHalfSizeMultiAtomonlyNewtonTri,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_MULTI_H
#define LMP_NPAIR_MULTI_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE, int ATOMONLY>
class NPairMulti : public NPair {
 public:
  NPairMulti(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
