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

#include "my_pool_chunk.h"

#include <cstdlib>

#if defined(LMP_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif

using namespace LAMMPS_NS;

/** \class LAMMPS_NS::MyPoolChunk
 *  \brief Templated class for storing chunks of datums in pages
 *
 * The size of the chunk may vary from call to call between the
 * *minchunk* and *maxchunk* setting.  Chunks may be returned
 * to the pool for re-use.  Chunks can be reserved in *nbin*
 * different sizes between *minchunk* and *maxchunk*.
 * The *chunksperpage* setting specifies how many chunks are stored
 * on any page and the *pagedelta* setting determines how many
 * pages are allocated in one go.  Pages are never freed, so they
 * can be re-used without re-allocation.
 *
 * \note
 * This is a template class with explicit instantiation. If the class
 * is used with a new data type a new explicit instantiation may need
 * to be added at the end of the file ``src/my_pool_chunk.cpp`` to
 * avoid symbol lookup errors. */

/** Create a class instance and set memory pool parameters
 *
 * \param  user_minchunk      Minimal chunk size
 * \param  user_maxchunk      Maximal chunk size
 * \param  user_nbin          Number of bins of different chunk sizes
 * \param  user_chunkperpage  Number of chunks per page
 * \param  user_pagedelta     Number of pages to allocate in one go */

template <class T>
MyPoolChunk<T>::MyPoolChunk(int user_minchunk, int user_maxchunk, int user_nbin,
                            int user_chunkperpage, int user_pagedelta)
{
  minchunk = user_minchunk;
  maxchunk = user_maxchunk;
  nbin = user_nbin;
  chunkperpage = user_chunkperpage;
  pagedelta = user_pagedelta;

  errorflag = 0;
  if (minchunk <= 0 || minchunk > maxchunk) errorflag = 1;
  if (user_nbin <= 0 || chunkperpage <= 0 || pagedelta <= 0) errorflag = 1;

  freehead = new int[nbin];
  chunksize = new int[nbin];
  if (!freehead || !chunksize) errorflag = 1;
  if (errorflag) return;

  // insure nbin*binsize spans minchunk to maxchunk inclusive

  binsize = (maxchunk - minchunk + 1) / nbin;
  if (minchunk + nbin * binsize <= maxchunk) binsize++;

  freelist = nullptr;
  for (int ibin = 0; ibin < nbin; ibin++) {
    freehead[ibin] = -1;
    chunksize[ibin] = minchunk + (ibin + 1) * binsize - 1;
    if (chunksize[ibin] > maxchunk) chunksize[ibin] = maxchunk;
  }

  ndatum = nchunk = 0;
  pages = nullptr;
  whichbin = nullptr;
  npage = 0;
}

/** Destroy class instance and free all allocated memory */
template <class T> MyPoolChunk<T>::~MyPoolChunk()
{
  delete[] freehead;
  delete[] chunksize;
  if (npage) {
    free(freelist);
    for (int i = 0; i < npage; i++) free(pages[i]);
    free(pages);
    free(whichbin);
  }
}

/** Return pointer/index of unused chunk of size maxchunk
 *
 * \param  index  Index of chunk in memory pool
 * \return        Pointer to requested chunk of storage */

template <class T> T *MyPoolChunk<T>::get(int &index)
{
  int ibin = nbin - 1;
  if (freehead[ibin] < 0) {
    allocate(ibin);
    if (errorflag) {
      index = -1;
      return nullptr;
    }
  }

  ndatum += maxchunk;
  nchunk++;
  index = freehead[ibin];
  int ipage = index / chunkperpage;
  int ientry = index % chunkperpage;
  freehead[ibin] = freelist[index];
  return &pages[ipage][ientry * chunksize[ibin]];
}

/** Return pointer/index of unused chunk of size N
 *
 * \param  n      Size of chunk
 * \param  index  Index of chunk in memory pool
 * \return        Pointer to requested chunk of storage */

template <class T> T *MyPoolChunk<T>::get(int n, int &index)
{
  if (n < minchunk || n > maxchunk) {
    errorflag = 3;
    index = -1;
    return nullptr;
  }

  int ibin = (n - minchunk) / binsize;
  if (freehead[ibin] < 0) {
    allocate(ibin);
    if (errorflag) {
      index = -1;
      return nullptr;
    }
  }

  ndatum += n;
  nchunk++;
  index = freehead[ibin];
  int ipage = index / chunkperpage;
  int ientry = index % chunkperpage;
  freehead[ibin] = freelist[index];
  return &pages[ipage][ientry * chunksize[ibin]];
}

/** Put indexed chunk back into memory pool via free list
 *
 * \param index  Memory chunk index returned by call to get() */

template <class T> void MyPoolChunk<T>::put(int index)
{
  if (index < 0) return;
  int ipage = index / chunkperpage;
  int ibin = whichbin[ipage];
  nchunk--;
  ndatum -= chunksize[ibin];
  freelist[index] = freehead[ibin];
  freehead[ibin] = index;
}

template <class T> void MyPoolChunk<T>::allocate(int ibin)
{
  int oldpage = npage;
  npage += pagedelta;
  freelist = (int *) realloc(freelist, sizeof(int) * npage * chunkperpage);
  pages = (T **) realloc(pages, sizeof(T *) * npage);
  whichbin = (int *) realloc(whichbin, sizeof(int) * npage);
  if (!freelist || !pages) {
    errorflag = 2;
    return;
  }

  // allocate pages with appropriate chunksize for ibin

  for (int i = oldpage; i < npage; i++) {
    whichbin[i] = ibin;
#if defined(LAMMPS_MEMALIGN)
    void *ptr;
    if (posix_memalign(&ptr, LAMMPS_MEMALIGN, sizeof(T) * chunkperpage * chunksize[ibin]))
      errorflag = 2;
    pages[i] = (T *) ptr;
#else
    pages[i] = (T *) malloc(sizeof(T) * chunkperpage * chunksize[ibin]);
    if (!pages[i]) errorflag = 2;
#endif
  }

  // reset free list for unused chunks on new pages

  freehead[ibin] = oldpage * chunkperpage;
  for (int i = freehead[ibin]; i < npage * chunkperpage; i++) freelist[i] = i + 1;
  freelist[npage * chunkperpage - 1] = -1;
}

/** Return total size of allocated pages
 *
 * \return total storage used in bytes */

template <class T> double MyPoolChunk<T>::size() const
{
  double bytes = (double) npage * chunkperpage * sizeof(int);
  bytes += (double) npage * sizeof(T *);
  bytes += (double) npage * sizeof(int);
  for (int i = 0; i < npage; ++i) bytes += (double) chunkperpage * chunksize[i] * sizeof(T);

  return bytes;
}

// explicit instantiations

namespace LAMMPS_NS {
template class MyPoolChunk<int>;
template class MyPoolChunk<double>;
}    // namespace LAMMPS_NS
