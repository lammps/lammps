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
typedef NPairNsq<0, 1, 0, 0> NPairFullNsq;
NPairStyle(full/nsq,
           NPairFullNsq,
           NP_FULL | NP_NSQ | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsq<1, 0, 0, 0> NPairHalfNsqNewtoff;
NPairStyle(half/nsq/newtoff,
           NPairHalfNsqNewtoff,
           NP_HALF | NP_NSQ | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsq<1, 1, 0, 0> NPairHalfNsqNewton;
NPairStyle(half/nsq/newton,
           NPairHalfNsqNewton,
           NP_HALF | NP_NSQ | NP_NEWTON | NP_ORTHO);

typedef NPairNsq<1, 1, 1, 0> NPairHalfNsqNewtonTri;
NPairStyle(half/nsq/newton/tri,
           NPairHalfNsqNewtonTri,
           NP_HALF | NP_NSQ | NP_NEWTON | NP_TRI);

typedef NPairNsq<0, 1, 0, 1> NPairFullSizeNsq;
NPairStyle(full/size/nsq,
           NPairFullSizeNsq,
           NP_FULL | NP_SIZE | NP_NSQ | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsq<1, 0, 0, 1> NPairHalfSizeNsqNewtoff;
NPairStyle(half/size/nsq/newtoff,
           NPairHalfSizeNsqNewtoff,
           NP_HALF | NP_SIZE | NP_NSQ | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsq<1, 1, 0, 1> NPairHalfSizeNsqNewton;
NPairStyle(half/size/nsq/newton,
           NPairHalfSizeNsqNewton,
           NP_HALF | NP_SIZE | NP_NSQ | NP_NEWTON | NP_ORTHO);

typedef NPairNsq<1, 1, 1, 1> NPairHalfSizeNsqNewtonTri;
NPairStyle(half/size/nsq/newton/tri,
           NPairHalfSizeNsqNewtonTri,
           NP_HALF | NP_SIZE | NP_NSQ | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_NSQ_H
#define LMP_NPAIR_NSQ_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE>
class NPairNsq : public NPair {
 public:
  NPairNsq(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
