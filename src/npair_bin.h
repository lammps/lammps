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
typedef NPairBin<0, 1, 0, 0, 0> NPairFullBin;
NPairStyle(full/bin,
           NPairFullBin,
           NP_FULL | NP_BIN | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 0, 0, 0, 0> NPairHalfBinNewtoff;
NPairStyle(half/bin/newtoff,
           NPairHalfBinNewtoff,
           NP_HALF | NP_BIN | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 1, 0, 0, 0> NPairHalfBinNewton;
NPairStyle(half/bin/newton,
           NPairHalfBinNewton,
           NP_HALF | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBin<1, 1, 1, 0, 0> NPairHalfBinNewtonTri;
NPairStyle(half/bin/newton/tri,
           NPairHalfBinNewtonTri,
           NP_HALF | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_TRI);

typedef NPairBin<0, 1, 0, 1, 0> NPairFullSizeBin;
NPairStyle(full/size/bin,
           NPairFullSizeBin,
           NP_FULL | NP_SIZE | NP_BIN | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 0, 0, 1, 0> NPairHalfSizeBinNewtoff;
NPairStyle(half/size/bin/newtoff,
           NPairHalfSizeBinNewtoff,
           NP_HALF | NP_SIZE | NP_BIN | NP_MOLONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 1, 0, 1, 0> NPairHalfSizeBinNewton;
NPairStyle(half/size/bin/newton,
           NPairHalfSizeBinNewton,
           NP_HALF | NP_SIZE | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBin<1, 1, 1, 1, 0> NPairHalfSizeBinNewtonTri;
NPairStyle(half/size/bin/newton/tri,
           NPairHalfSizeBinNewtonTri,
           NP_HALF | NP_SIZE | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_TRI);

typedef NPairBin<0, 1, 0, 0, 1> NPairFullBinAtomonly;
NPairStyle(full/bin/atomonly,
           NPairFullBinAtomonly,
           NP_FULL | NP_BIN | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 0, 0, 0, 1> NPairHalfBinAtomonlyNewtoff;
NPairStyle(half/bin/atomonly/newtoff,
           NPairHalfBinAtomonlyNewtoff,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 1, 0, 0, 1> NPairHalfBinAtomonlyNewton;
NPairStyle(half/bin/atomonly/newton,
           NPairHalfBinAtomonlyNewton,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBin<1, 1, 1, 0, 1> NPairHalfBinAtomonlyNewtonTri;
NPairStyle(half/bin/atomonly/newton/tri,
           NPairHalfBinAtomonlyNewtonTri,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_TRI);

typedef NPairBin<0, 1, 0, 1, 1> NPairFullSizeBinAtomonly;
NPairStyle(full/size/bin/atomonly,
           NPairFullSizeBinAtomonly,
           NP_FULL | NP_SIZE | NP_BIN | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 0, 0, 1, 1> NPairHalfSizeBinAtomonlyNewtoff;
NPairStyle(half/size/bin/atomonly/newtoff,
           NPairHalfSizeBinAtomonlyNewtoff,
           NP_HALF | NP_SIZE | NP_BIN | NP_ATOMONLY | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBin<1, 1, 0, 1, 1> NPairHalfSizeBinAtomonlyNewton;
NPairStyle(half/size/bin/atomonly/newton,
           NPairHalfSizeBinAtomonlyNewton,
           NP_HALF | NP_SIZE | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);

typedef NPairBin<1, 1, 1, 1, 1> NPairHalfSizeBinAtomonlyNewtonTri;
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
