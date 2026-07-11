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
using NPairFullBin = NPairBin<0, 1, 0, 0, 0>;
NPairStyle(full/bin,
           NPairFullBin,
           NP_FULL | NP_BIN | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinNewtoff = NPairBin<1, 0, 0, 0, 0>;
NPairStyle(half/bin/newtoff,
           NPairHalfBinNewtoff,
           NP_HALF | NP_BIN | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinNewton = NPairBin<1, 1, 0, 0, 0>;
NPairStyle(half/bin/newton,
           NPairHalfBinNewton,
           NP_HALF | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfBinNewtonTri = NPairBin<1, 1, 1, 0, 0>;
NPairStyle(half/bin/newton/tri,
           NPairHalfBinNewtonTri,
           NP_HALF | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_TRI);

using NPairFullSizeBin = NPairBin<0, 1, 0, 1, 0>;
NPairStyle(full/size/bin,
           NPairFullSizeBin,
           NP_FULL | NP_SIZE | NP_BIN | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinNewtoff = NPairBin<1, 0, 0, 1, 0>;
NPairStyle(half/size/bin/newtoff,
           NPairHalfSizeBinNewtoff,
           NP_HALF | NP_SIZE | NP_BIN | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinNewton = NPairBin<1, 1, 0, 1, 0>;
NPairStyle(half/size/bin/newton,
           NPairHalfSizeBinNewton,
           NP_HALF | NP_SIZE | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeBinNewtonTri = NPairBin<1, 1, 1, 1, 0>;
NPairStyle(half/size/bin/newton/tri,
           NPairHalfSizeBinNewtonTri,
           NP_HALF | NP_SIZE | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_TRI);

using NPairFullBinAtomonly = NPairBin<0, 1, 0, 0, 1>;
NPairStyle(full/bin/atomonly,
           NPairFullBinAtomonly,
           NP_FULL | NP_BIN | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinAtomonlyNewtoff = NPairBin<1, 0, 0, 0, 1>;
NPairStyle(half/bin/atomonly/newtoff,
           NPairHalfBinAtomonlyNewtoff,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfBinAtomonlyNewton = NPairBin<1, 1, 0, 0, 1>;
NPairStyle(half/bin/atomonly/newton,
           NPairHalfBinAtomonlyNewton,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfBinAtomonlyNewtonTri = NPairBin<1, 1, 1, 0, 1>;
NPairStyle(half/bin/atomonly/newton/tri,
           NPairHalfBinAtomonlyNewtonTri,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_TRI);

using NPairFullSizeBinAtomonly = NPairBin<0, 1, 0, 1, 1>;
NPairStyle(full/size/bin/atomonly,
           NPairFullSizeBinAtomonly,
           NP_FULL | NP_SIZE | NP_BIN | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinAtomonlyNewtoff = NPairBin<1, 0, 0, 1, 1>;
NPairStyle(half/size/bin/atomonly/newtoff,
           NPairHalfSizeBinAtomonlyNewtoff,
           NP_HALF | NP_SIZE | NP_BIN | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

using NPairHalfSizeBinAtomonlyNewton = NPairBin<1, 1, 0, 1, 1>;
NPairStyle(half/size/bin/atomonly/newton,
           NPairHalfSizeBinAtomonlyNewton,
           NP_HALF | NP_SIZE | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

using NPairHalfSizeBinAtomonlyNewtonTri = NPairBin<1, 1, 1, 1, 1>;
NPairStyle(half/size/bin/atomonly/newton/tri,
           NPairHalfSizeBinAtomonlyNewtonTri,
           NP_HALF | NP_SIZE | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_BIN_H
#define LMP_NPAIR_BIN_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE, int ATOMONLY>
class NPairBin : public NPair {
 public:
  NPairBin(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
