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

#include "my_page.h"

#if defined(LMP_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif

using namespace LAMMPS_NS;

/** \class LAMMPS_NS::MyPage
 * \brief Templated class for storing chunks of datums in pages.
 *
 * The size of the chunk may vary from call to call, but must be
 * less or equal than the *maxchunk* setting.
 * The chunks are not returnable like with malloc() (i.e. you cannot
 * call free() on them individually).  One can only reset and start over.
 * The purpose of this class is to replace many small memory allocations
 * via malloc() with a few large ones.  Since the pages are never freed
 * until the class is re-initialized, they can be re-used without having
 * to re-allocate them by calling the reset() method.
 *
 * The settings *maxchunk*, *pagesize*, and *pagedelta* control
 * the memory allocation strategy.  The *maxchunk* value represents
 * the expected largest number of items per chunk.  If there is
 * less space left on the current page, a new page is allocated
 * for the next chunk.  The *pagesize* value represents how many
 * items can fit on a single page.  It should have space for multiple
 * chunks of size *maxchunk*.  The combination of these two
 * parameters determines how much memory is wasted by either switching
 * to the next page too soon or allocating too large pages that never
 * get properly used.  It is an error, if a requested chunk is larger
 * than *maxchunk*.  The *pagedelta* parameter determines how many
 * pages are allocated in one go.  In combination with the *pagesize*
 * setting, this determines how often blocks of memory get allocated
 * (fewer allocations will result in faster execution).
 *
 * \note
 * This is a template class with explicit instantiation. If the class
 * is used with a new data type a new explicit instantiation may need to
 * be added at the end of the file ``src/my_page.cpp`` to avoid symbol
 * lookup errors. */

/** Create a class instance
 *
 *  Need to call init() before use to define allocation settings */

template <class T>
MyPage<T>::MyPage() :
    ndatum(0), nchunk(0), pages(nullptr), page(nullptr), npage(0), ipage(-1), index(-1),
    maxchunk(-1), pagesize(-1), pagedelta(1), errorflag(0){};

template <class T> MyPage<T>::~MyPage()
{
  deallocate();
}

/** (Re-)initialize the set of pages and allocation parameters.
 *
 * This also frees all previously allocated storage and allocates
 * the first page(s).
 *
 * \param  user_maxchunk   Expected maximum number of items for one chunk
 * \param  user_pagesize   Number of items on a single memory page
 * \param  user_pagedelta  Number of pages to allocate with one malloc
 * \return                 1 if there were invalid parameters, 2 if there was an allocation error or 0 if successful */

template <class T> int MyPage<T>::init(int user_maxchunk, int user_pagesize, int user_pagedelta)
{
  maxchunk = user_maxchunk;
  pagesize = user_pagesize;
  pagedelta = user_pagedelta;

  if (maxchunk <= 0 || pagesize <= 0 || pagedelta <= 0) return 1;
  if (maxchunk > pagesize) return 1;

  // free storage if re-initialized

  deallocate();

  // initial page allocation

  allocate();
  if (errorflag) return 2;
  reset();
  return 0;
}

/** Pointer to location that can store N items.
 *
 * This will allocate more pages as needed.
 * If the parameter *N* is larger than the *maxchunk*
 * setting an error is flagged.
 *
 * \param  n  number of items for which storage is requested
 * \return    memory location or null pointer, if error or allocation failed */

template <class T> T *MyPage<T>::get(int n)
{
  if (n > maxchunk) {
    errorflag = 1;
    return nullptr;
  }
  ndatum += n;
  nchunk++;

  // return pointer from current page
  if (index + n <= pagesize) {
    int start = index;
    index += n;
    return &page[start];
  }

  // allocate new page
  ipage++;
  if (ipage == npage) {
    allocate();
    if (errorflag) return nullptr;
  }
  page = pages[ipage];
  index = n;
  return &page[0];
}

/** Reset state of memory pool without freeing any memory */

template <class T> void MyPage<T>::reset()
{
  ndatum = nchunk = 0;
  index = ipage = 0;
  page = (pages != nullptr) ? pages[ipage] : nullptr;
  errorflag = 0;
}

/* ---------------------------------------------------------------------- */

template <class T> void MyPage<T>::allocate()
{
  npage += pagedelta;
  pages = (T **) realloc(pages, npage * sizeof(T *));
  if (!pages) {
    errorflag = 2;
    return;
  }

  for (int i = npage - pagedelta; i < npage; i++) {
#if defined(LAMMPS_MEMALIGN)
    void *ptr;
    if (posix_memalign(&ptr, LAMMPS_MEMALIGN, pagesize * sizeof(T))) errorflag = 2;
    pages[i] = (T *) ptr;
#else
    pages[i] = (T *) malloc(pagesize * sizeof(T));
    if (!pages[i]) errorflag = 2;
#endif
  }
}

/** Free all allocated pages of this class instance */

template <class T> void MyPage<T>::deallocate()
{
  reset();
  for (int i = 0; i < npage; i++) free(pages[i]);
  free(pages);
  pages = nullptr;
  npage = 0;
}

// explicit instantiations

namespace LAMMPS_NS {
template class MyPage<int>;
template class MyPage<long>;
template class MyPage<long long>;
template class MyPage<double>;
template class MyPage<HyperOneCoeff>;
}    // namespace LAMMPS_NS
