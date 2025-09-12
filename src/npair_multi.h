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
using NPairFullMulti = NPairMulti<0, 1, 0, 0, 0>;
NPairStyle(full/multi,
           NPairFullMulti,
           NP_FULL | NP_MULTI | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiNewtoff = NPairMulti<1, 0, 0, 0, 0>;
NPairStyle(half/multi/newtoff,
           NPairHalfMultiNewtoff,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiNewton = NPairMulti<1, 1, 0, 0, 0>;
NPairStyle(half/multi/newton,
           NPairHalfMultiNewton,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfMultiNewtonTri = NPairMulti<1, 1, 1, 0, 0>;
NPairStyle(half/multi/newton/tri,
           NPairHalfMultiNewtonTri,
           NP_HALF | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_TRI);

using NPairFullSizeMulti = NPairMulti<0, 1, 0, 1, 0>;
NPairStyle(full/size/multi,
           NPairFullSizeMulti,
           NP_FULL | NP_SIZE | NP_MULTI | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiNewtoff = NPairMulti<1, 0, 0, 1, 0>;
NPairStyle(half/size/multi/newtoff,
           NPairHalfSizeMultiNewtoff,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiNewton = NPairMulti<1, 1, 0, 1, 0>;
NPairStyle(half/size/multi/newton,
           NPairHalfSizeMultiNewton,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeMultiNewtonTri = NPairMulti<1, 1, 1, 1, 0>;
NPairStyle(half/size/multi/newton/tri,
           NPairHalfSizeMultiNewtonTri,
           NP_HALF | NP_SIZE | NP_MULTI | NP_MOLONLY | NP_NEWTON | NP_TRI);

using NPairFullMultiAtomonly = NPairMulti<0, 1, 0, 0, 1>;
NPairStyle(full/multi/atomonly,
           NPairFullMultiAtomonly,
           NP_FULL | NP_MULTI | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiAtomonlyNewtoff = NPairMulti<1, 0, 0, 0, 1>;
NPairStyle(half/multi/atomonly/newtoff,
           NPairHalfMultiAtomonlyNewtoff,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfMultiAtomonlyNewton = NPairMulti<1, 1, 0, 0, 1>;
NPairStyle(half/multi/atomonly/newton,
           NPairHalfMultiAtomonlyNewton,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfMultiAtomonlyNewtonTri = NPairMulti<1, 1, 1, 0, 1>;
NPairStyle(half/multi/atomonly/newton/tri,
           NPairHalfMultiAtomonlyNewtonTri,
           NP_HALF | NP_MULTI | NP_ATOMONLY | NP_NEWTON | NP_TRI);

using NPairFullSizeMultiAtomonly = NPairMulti<0, 1, 0, 1, 1>;
NPairStyle(full/size/multi/atomonly,
           NPairFullSizeMultiAtomonly,
           NP_FULL | NP_SIZE | NP_MULTI | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiAtomonlyNewtoff = NPairMulti<1, 0, 0, 1, 1>;
NPairStyle(half/size/multi/atomonly/newtoff,
           NPairHalfSizeMultiAtomonlyNewtoff,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeMultiAtomonlyNewton = NPairMulti<1, 1, 0, 1, 1>;
NPairStyle(half/size/multi/atomonly/newton,
           NPairHalfSizeMultiAtomonlyNewton,
           NP_HALF | NP_SIZE | NP_MULTI | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeMultiAtomonlyNewtonTri = NPairMulti<1, 1, 1, 1, 1>;
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
