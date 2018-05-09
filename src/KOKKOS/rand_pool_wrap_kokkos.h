/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef RAND_POOL_WRAP_H
#define RAND_POOL_WRAP_H

#include "pointers.h"
#include "kokkos_type.h"
#include "random_mars.h"
#include "error.h"

namespace LAMMPS_NS {

struct RandWrap {
  class RanMars* rng;

  KOKKOS_INLINE_FUNCTION
  RandWrap() {
    rng = NULL;
  }

  KOKKOS_INLINE_FUNCTION
  double drand() {
    return rng->uniform();
  }

  KOKKOS_INLINE_FUNCTION
  double normal() {
    return rng->gaussian();
  }
};

class RandPoolWrap : protected Pointers {
 public:
  RandPoolWrap(int, class LAMMPS *);
  ~RandPoolWrap();
  void destroy();
  void init(RanMars*, int);

  KOKKOS_INLINE_FUNCTION
  RandWrap get_state() const
  {
#ifdef KOKKOS_HAVE_CUDA
    error->all(FLERR,"Cannot use Marsaglia RNG with GPUs");
#endif

    RandWrap rand_wrap;

    typedef Kokkos::Experimental::UniqueToken<
      LMPHostType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;

    unique_token_type unique_token;
    int tid = (int) unique_token.acquire();
    rand_wrap.rng = random_thr[tid];
    unique_token.release(tid);

    return rand_wrap;
  }

  KOKKOS_INLINE_FUNCTION
  void free_state(RandWrap) const
  {

  }

 private:
  class RanMars **random_thr;
  int nthreads;
};

}

#endif

/* ERROR/WARNING messages:

*/
