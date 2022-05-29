// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NPAIR_INTEL_H
#define LMP_NPAIR_INTEL_H

#include "npair.h"
#include "domain.h"
#include "fix_intel.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef LMP_USE_AVXCD
#include "intel_simd.h"
#endif

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push,target(mic))
#endif

#define ofind_special(which, special, nspecial, i, tag)               \
{                                                                     \
  which = 0;                                                          \
  const int n1 = nspecial[i * 3];                                     \
  const int n2 = nspecial[i * 3 + 1];                                 \
  const int n3 = nspecial[i * 3 + 2];                                 \
  const tagint *sptr = special + i * maxspecial;                      \
  for (int s = 0; s < n3; s++) {                                      \
    if (sptr[s] == tag) {                                             \
      if (s < n1) {                                                   \
        which = 1;                                                    \
      } else if (s < n2) {                                            \
        which = 2;                                                    \
      } else {                                                        \
        which = 3;                                                    \
      }                                                               \
    }                                                                 \
  }                                                                   \
}

#define ominimum_image_check(answer, dx, dy, dz)                      \
{                                                                     \
  answer = 0;                                                         \
  if (xperiodic && fabs(dx) > xprd_half) answer = 1;                  \
  if (yperiodic && fabs(dy) > yprd_half) answer = 1;                  \
  if (zperiodic && fabs(dz) > zprd_half) answer = 1;                  \
}

#define dminimum_image_check(answer, dx, dy, dz)                      \
{                                                                     \
  answer = 0;                                                         \
  if (domain->xperiodic && fabs(dx) > domain->xprd_half) answer = 1;  \
  if (domain->yperiodic && fabs(dy) > domain->yprd_half) answer = 1;  \
  if (domain->zperiodic && fabs(dz) > domain->zprd_half) answer = 1;  \
}

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

namespace LAMMPS_NS {

class NPairIntel : public NPair {
 public:
  NPairIntel(class LAMMPS *);
  #ifdef _LMP_INTEL_OFFLOAD
  ~NPairIntel() override;
  #endif
  void copy_neighbor_info() override;

  #ifdef _LMP_INTEL_OFFLOAD
  void grow_stencil();
  #endif

 protected:
  FixIntel *_fix;

  template <class flt_t, class acc_t>
  void copy_cutsq_info(IntelBuffers<flt_t,acc_t> *);

  template <class flt_t, class acc_t, int, int, int, int, int>
  void bin_newton(const int, NeighList *, IntelBuffers<flt_t,acc_t> *,
                  const int, const int, const int offload_end = 0);

  #ifdef _LMP_INTEL_OFFLOAD
  int _cop;
  int *_off_map_stencil;
  #endif
};

}

#endif


