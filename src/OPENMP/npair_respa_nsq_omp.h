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
typedef NPairRespaNsqOmp<0,0> NPairHalfRespaNsqNewtoffOmp;
NPairStyle(half/respa/nsq/newtoff/omp,
           NPairHalfRespaNsqNewtoffOmp,
           NP_HALF | NP_RESPA | NP_NSQ | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairRespaNsqOmp<1,0> NPairHalfRespaNsqNewtonOmp;
NPairStyle(half/respa/nsq/newton/omp,
           NPairHalfRespaNsqNewtonOmp,
           NP_HALF | NP_RESPA | NP_NSQ | NP_OMP | NP_NEWTON | NP_ORTHO);

typedef NPairRespaNsqOmp<1,1> NPairHalfRespaNsqNewtonTriOmp;
NPairStyle(half/respa/nsq/newton/tri/omp,
           NPairHalfRespaNsqNewtonTriOmp,
           NP_HALF | NP_RESPA | NP_NSQ | NP_OMP | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_RESPA_NSQ_OMP_H
#define LMP_NPAIR_RESPA_NSQ_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int NEWTON, int TRI>
class NPairRespaNsqOmp : public NPair {
 public:
  NPairRespaNsqOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
