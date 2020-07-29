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

#include "memory.h"
#include <cstdlib>
#include "error.h"
#include "fmt/format.h"

#if defined(LMP_USER_INTEL) && defined(__INTEL_COMPILER)
#ifndef LMP_INTEL_NO_TBB
#define LMP_USE_TBB_ALLOCATOR
#include "tbb/scalable_allocator.h"
#else
#include <cstring>
#if defined(__APPLE__)
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#endif
#endif

#if defined(LMP_USER_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Memory::Memory(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   safe malloc
------------------------------------------------------------------------- */

void *Memory::smalloc(bigint nbytes, const char *name)
{
  if (nbytes == 0) return NULL;

#if defined(LAMMPS_MEMALIGN)
  void *ptr;

#if defined(LMP_USE_TBB_ALLOCATOR)
  ptr = scalable_aligned_malloc(nbytes, LAMMPS_MEMALIGN);
#else
  int retval = posix_memalign(&ptr, LAMMPS_MEMALIGN, nbytes);
  if (retval) ptr = NULL;
#endif

#else
  void *ptr = malloc(nbytes);
#endif
  if (ptr == NULL)
    error->one(FLERR,fmt::format("Failed to allocate {} bytes for array {}",
                                 nbytes,name));
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
  if (nbytes == 0) {
    destroy(ptr);
    return NULL;
  }

#if defined(LMP_USE_TBB_ALLOCATOR)
  ptr = scalable_aligned_realloc(ptr, nbytes, LAMMPS_MEMALIGN);
#elif defined(LMP_INTEL_NO_TBB) && defined(LAMMPS_MEMALIGN) && \
      defined(__INTEL_COMPILER)

  ptr = realloc(ptr, nbytes);
  uintptr_t offset = ((uintptr_t)(const void *)(ptr)) % LAMMPS_MEMALIGN;
  if (offset) {
    void *optr = ptr;
    ptr = smalloc(nbytes, name);
#if defined(__APPLE__)
    memcpy(ptr, optr, MIN(nbytes,malloc_size(optr)));
#elif defined(_WIN32) || defined(__MINGW32__)
    memcpy(ptr, optr, MIN(nbytes,_msize(optr)));
#else
    memcpy(ptr, optr, MIN(nbytes,malloc_usable_size(optr)));
#endif
    free(optr);
  }
#else
  ptr = realloc(ptr,nbytes);
#endif
  if (ptr == NULL)
    error->one(FLERR,fmt::format("Failed to reallocate {} bytes for array {}",
                                 nbytes,name));
  return ptr;
}

/* ----------------------------------------------------------------------
   safe free
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  #if defined(LMP_USE_TBB_ALLOCATOR)
  scalable_aligned_free(ptr);
  #else
  free(ptr);
  #endif
}

/* ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(const char *name)
{
  error->one(FLERR,fmt::format("Cannot create/grow a vector/array of "
                               "pointers for {}",name));
}
