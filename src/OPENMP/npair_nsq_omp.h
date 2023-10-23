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

typedef NPairNsqOmp<0, 1, 0, 0> NPairFullNsqOmp;
NPairStyle(full/nsq/omp,
           NPairFullNsqOmp,
           NP_FULL | NP_NSQ | NP_OMP | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsqOmp<1, 0, 0, 0> NPairHalfNsqNewtoffOmp;
NPairStyle(half/nsq/newtoff/omp,
           NPairHalfNsqNewtoffOmp,
           NP_HALF | NP_NSQ | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsqOmp<1, 1, 0, 0> NPairHalfNsqNewtonOmp;
NPairStyle(half/nsq/newton/omp,
           NPairHalfNsqNewtonOmp,
           NP_HALF | NP_NSQ | NP_OMP | NP_NEWTON | NP_ORTHO);

typedef NPairNsqOmp<1, 1, 1, 0> NPairHalfNsqNewtonTriOmp;
NPairStyle(half/nsq/newton/tri/omp,
           NPairHalfNsqNewtonTriOmp,
           NP_HALF | NP_NSQ | NP_OMP | NP_NEWTON | NP_TRI);

typedef NPairNsqOmp<0, 1, 0, 1> NPairFullSizeNsqOmp;
NPairStyle(full/size/nsq/omp,
           NPairFullSizeNsqOmp,
           NP_FULL | NP_SIZE | NP_NSQ | NP_OMP | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsqOmp<1, 0, 0, 1> NPairHalfSizeNsqNewtoffOmp;
NPairStyle(half/size/nsq/newtoff/omp,
           NPairHalfSizeNsqNewtoffOmp,
           NP_HALF | NP_SIZE | NP_NSQ | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairNsqOmp<1, 1, 0, 1> NPairHalfSizeNsqNewtonOmp;
NPairStyle(half/size/nsq/newton/omp,
           NPairHalfSizeNsqNewtonOmp,
           NP_HALF | NP_SIZE | NP_NSQ | NP_OMP | NP_NEWTON | NP_ORTHO);

typedef NPairNsqOmp<1, 1, 1, 1> NPairHalfSizeNsqNewtonTriOmp;
NPairStyle(half/size/nsq/newton/tri/omp,
           NPairHalfSizeNsqNewtonTriOmp,
           NP_HALF | NP_SIZE | NP_NSQ | NP_OMP | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_NSQ_OMP_H
#define LMP_NPAIR_NSQ_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE>
class NPairNsqOmp : public NPair {
 public:
  NPairNsqOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
