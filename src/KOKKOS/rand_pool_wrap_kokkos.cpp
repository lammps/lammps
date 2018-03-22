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

#include "comm.h"
#include "rand_pool_wrap_kokkos.h"
#include "lammps.h"
#include "kokkos.h"
#include "random_mars.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RandPoolWrap::RandPoolWrap(int, LAMMPS *lmp) : Pointers(lmp)
{
  random_thr =  NULL;
  nthreads = lmp->kokkos->num_threads;
}

/* ---------------------------------------------------------------------- */

RandPoolWrap::~RandPoolWrap()
{

}

void RandPoolWrap::destroy()
{
  if (random_thr) {
    for (int i=1; i < nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
    random_thr = NULL;
  }
}

void RandPoolWrap::init(RanMars* random, int seed)
{
  // deallocate pool of RNGs
  if (random_thr) {
    for (int i=1; i < this->nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
  }

  // allocate pool of RNGs
  // generate a random number generator instance for
  // all threads != 0. make sure we use unique seeds.
  nthreads = lmp->kokkos->num_threads;
  random_thr = new RanMars*[nthreads];
  for (int tid = 1; tid < nthreads; ++tid) {
    random_thr[tid] = new RanMars(lmp, seed + comm->me
                                  + comm->nprocs*tid);
  }

  // to ensure full compatibility with the serial style
  // we use the serial random number generator instance for thread 0
  random_thr[0] = random;
}
