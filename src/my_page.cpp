/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
usage:
  request one datum at a time, repeat, clear
  request chunks of datums in each get() or vget(), repeat, clear
  chunk size can vary from request to request
  chunk size can be known in advance or registered after usage via vgot()
inputs:
   template T = one datum, e.g. int, double, struct, int[3]
     for int[3], access datum as ivec[i][2]
methods:
   T *get() = return ptr to one datum
   T *get(N) = return ptr to N datums, N < maxchunk required
   T *vget() = return ptr to maxchunk datums, use as needed, then call vgot()
     all gets return NULL if error encountered
   vgot(N) = used N datums of previous vget(), N < maxchunk required
   void init(maxchunk, pagesize, pagedelta)
     define allocation params and allocate first page(s)
     call right after constructor
       can call again to reset allocation params and free previous pages
     maxchunk = max # of datums in one chunk, default = 1
     pagesize = # of datums in one page, default = 1024
       should be big enough to store multiple chunks
     pagedelta = # of pages to allocate at a time, default = 1
     return 1 if bad params
   void reset() = clear pages w/out freeing
   int size() = return total size of allocated pages in bytes
   int status() = return error status
     0 = ok, 1 = chunksize > maxchunk, 2 = allocation error
------------------------------------------------------------------------- */

#include "my_page.h"

#include <cstdlib>

#if defined(LMP_USER_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif

using namespace LAMMPS_NS;

/** \class LAMMPS_NS::MyPage
 * \brief Templated class for storing chunks of datums in pages.
 *
 * The chunks are not returnable (i.e. cannot be freed individually).
 * One can only reset and start over.  The purpose of this
 * class is to replace many small mallocs with a few large
 * mallocs. Since the pages are never freed, they can be re-used
 * without having to re-allocate them.
 *
 * The settings *maxchunk*, *pagesize*, and *pagedelta* contol
 * the memory allocation strategy. The *maxchunk* value represents
 * the expected largest number of items per chunk. If there is
 * less space left on the current page, a new page is allocated
 * for the next chunk.  The *pagesize* value represents how many
 * items can fit on a single page. It should have space for multiple
 * chunks of size *maxchunk*.  The combination of these two
 * parameters determines how much memory is wasted by either switching
 * to the next page too soon or allocating too large pages that never
 * get properly used.  It is an error, if a requested chunk is larger
 * than *maxchunk*.  The *pagedelta* parameter determines how many
 * pages are allocated in one go.  In combination with the *pagesize*
 * setting, this determines how often blocks of memory get allocated
 * (fewer allocations will result in faster execution).
 *
 * This class is the "workhorse" for the memory management of
 * neighbor lists. */

/** Create a class instance
 *
 *  Need to call init() before use to define allocation settings */

template <class T>
MyPage<T>::MyPage() : ndatum(0), nchunk(0), pages(nullptr), page(nullptr),
                      npage(0), ipage(-1), index(-1), maxchunk(-1),
                      pagesize(-1), pagedelta(1), errorflag(0) {};

/** Free all allocated pages of this class instance */

template <class T>
MyPage<T>::~MyPage() {
  for (int i = 0; i < npage; i++) free(pages[i]);
  free(pages);
}

/** (Re-)initialize the set of pages and allocation parameters.
 *
 * This also frees all previously allocated storage and allocates
 * the first page(s).
 *
 * \param  user_maxchunk    Expected maximum number of items for one chunk
 * \param  user_pagesize    Number of items on a single memory page
 * \param  user_page_delta  Number of pages to allocate with one malloc
 * \return                  1 if there was an error or 0 if successful */

template<class T>
int MyPage<T>::init(int user_maxchunk, int user_pagesize,
           int user_pagedelta) {
    maxchunk = user_maxchunk;
    pagesize = user_pagesize;
    pagedelta = user_pagedelta;

    if (maxchunk <= 0 || pagesize <= 0 || pagedelta <= 0) return 1;
    if (maxchunk > pagesize) return 1;

    // free any previously allocated pages

    for (int i = 0; i < npage; i++) free(pages[i]);
    free(pages);

    // initial page allocation

    ndatum = nchunk = 0;
    pages = NULL;
    npage = 0;
    allocate();
    if (errorflag) return 2;
    ipage = index = 0;
    page = pages[ipage];
    return 0;
  }

/** Pointer to location that can store one item.
 *
 * This will allocate more pages as needed.
 *
 * \return  memory location or null pointer, if memory allocation failed */

template <class T>
T *MyPage<T>::get() {
  ndatum++;
  nchunk++;
  if (index < pagesize) return &page[index++];
  ipage++;
  if (ipage == npage) {
    allocate();
    if (errorflag) return NULL;
  }
  page = pages[ipage];
  index = 0;
  return &page[index++];
}

/** Pointer to location that can store N items.
 *
 * This will allocate more pages as needed.
 * If the parameter *N* is larger than the *maxchunk*
 * setting an error is flagged.
 *
 * \param  n  number of items for which storage is requested
 * \return    memory location or null pointer, if error or allocation failed */

template <class T>
T *MyPage<T>::get(int n) {
  if (n > maxchunk) {
    errorflag = 1;
    return NULL;
  }
  ndatum += n;
  nchunk++;
  if (index+n <= pagesize) {
    int start = index;
    index += n;
    return &page[start];
  }
  ipage++;
  if (ipage == npage) {
    allocate();
    if (errorflag) return NULL;
  }
  page = pages[ipage];
  index = n;
  return &page[0];
}

/** Get pointer to location that can store *maxchunk* items.
 *
 * This will return the same pointer as the previous call to
 * this function unless vgot() is called afterwards to record
 * how many items of the chunk were actually used.
 *
 * \return pointer to chunk of memory or null pointer if run out of memory */

template <class T>
T *MyPage<T>::vget() {
  if (index+maxchunk <= pagesize) return &page[index];
  ipage++;
  if (ipage == npage) {
    allocate();
    if (errorflag) return NULL;
  }
  page = pages[ipage];
  index = 0;
  return &page[index];
}

/** Mark *N* items as used of the chunk reserved with a preceding call to vget().
 *
 * This will advance the internal pointer inside the current memory page.
 * It is not necessary to call this function for *N* = 0, that is the reserved
 * storage was not used.  A following call to vget() will then reserve the
 * same location again.  It is an error if *N* > *maxchunk*.
 *
 * \param  n  Number of iterms used in previously reserved chunk */

template <class T>
void MyPage<T>::vgot(int n) {
  if (n > maxchunk) errorflag = 1;
  ndatum += n;
  nchunk++;
  index += n;
}

/** Reset state of memory pool without freeing any memory */

template <class T>
void MyPage<T>::reset() {
  ndatum = nchunk = 0;
  index = ipage = 0;
  page = pages[ipage];
}

/* ---------------------------------------------------------------------- */

template <class T>
void MyPage<T>::allocate() {
  npage += pagedelta;
  pages = (T **) realloc(pages,npage*sizeof(T *));
  if (!pages) {
    errorflag = 2;
    return;
  }

  for (int i = npage-pagedelta; i < npage; i++) {
#if defined(LAMMPS_MEMALIGN)
    void *ptr;
    if (posix_memalign(&ptr, LAMMPS_MEMALIGN, pagesize*sizeof(T)))
      errorflag = 2;
    pages[i] = (T *) ptr;
#else
    pages[i] = (T *) malloc(pagesize*sizeof(T));
    if (!pages[i]) errorflag = 2;
#endif
  }
}

// explicit instantiations

namespace LAMMPS_NS {
  template class MyPage<int>;
  template class MyPage<long>;
  template class MyPage<long long>;
  template class MyPage<double>;
  template class MyPage<HyperOneCoeff>;
}
