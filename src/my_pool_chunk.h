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
MyPoolChunk = templated class for storing chunks of datums in pages
  chunks can be returned to pool for reuse
  chunks come in nbin different fixed sizes so can reuse
  replaces many small mallocs with a few large mallocs
  pages are never freed, so can reuse w/out reallocs
usage:
  continously get() and put() chunks as needed
  NOTE: could add a clear() if retain info on mapping of pages to bins
inputs:
   template T = one datum, e.g. int, double, struct
   minchunk = min # of datums in one chunk, def = 1
   maxchunk = max # of datums in one chunk, def = 1
   nbin = # of bins between minchunk and maxchunk
   chunkperpage = # of chunks in one page, def = 1024
   pagedelta = # of pages to allocate at a time, def = 1
methods:
   T *get(index) = return ptr/index to unused chunk of size maxchunk
   T *get(N,index) = return ptr/index to unused chunk of size N
                     minchunk <= N <= maxchunk required
   put(index) = return indexed chunk to pool (same index returned by get)
   int size() = return total size of allocated pages in bytes
public variables:
   ndatum = total # of stored datums
   nchunk = total # of stored chunks
   size = total size of all allocated pages in daums
   errorflag = flag for various error conditions
------------------------------------------------------------------------- */

#ifndef LAMMPS_MY_POOL_CHUNK_H
#define LAMMPS_MY_POOL_CHUNK_H

#if defined(LMP_USER_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif

#include <cstdlib>

namespace LAMMPS_NS {

template<class T>
class MyPoolChunk {
 public:
  int ndatum;      // total # of stored datums
  int nchunk;      // total # of stored chunks
  int size;        // total size of all allocated pages in datums
  int errorflag;   // flag > 1 if error has occurred
                   // 1 = invalid inputs
                   // 2 = memory allocation error
                   // 3 = chunk size exceeded maxchunk

  MyPoolChunk(int user_minchunk = 1, int user_maxchunk = 1, int user_nbin = 1,
              int user_chunkperpage = 1024, int user_pagedelta = 1) {
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

    binsize = (maxchunk-minchunk+1) / nbin;
    if (minchunk + nbin*binsize <= maxchunk) binsize++;

    freelist = NULL;
    for (int ibin = 0; ibin < nbin; ibin++) {
      freehead[ibin] = -1;
      chunksize[ibin] = minchunk + (ibin+1)*binsize - 1;
      if (chunksize[ibin] > maxchunk) chunksize[ibin] = maxchunk;
    }

    ndatum = nchunk = size = 0;
    pages = NULL;
    whichbin = NULL;
    npage = 0;
  }

  // free all allocated memory

  ~MyPoolChunk() {
    delete [] freehead;
    delete [] chunksize;
    if (npage) {
      free(freelist);
      for (int i = 0; i < npage; i++) free(pages[i]);
      free(pages);
      free(whichbin);
    }
  }

  // return pointer/index of unused chunk of size maxchunk

  T *get(int &index) {
    int ibin = nbin-1;
    if (freehead[ibin] < 0) {
      allocate(ibin);
      if (errorflag) return NULL;
    }

    ndatum += maxchunk;
    nchunk++;
    index = freehead[ibin];
    int ipage = index/chunkperpage;
    int ientry = index % chunkperpage;
    freehead[ibin] = freelist[index];
    return &pages[ipage][ientry*chunksize[ibin]];
  }

  // return pointer/index of unused chunk of size N

  T *get(int n, int &index) {
    if (n < minchunk || n > maxchunk) {
      errorflag = 3;
      return NULL;
    }

    int ibin = (n-minchunk) / binsize;
    if (freehead[ibin] < 0) {
      allocate(ibin);
      if (errorflag) return NULL;
    }

    ndatum += n;
    nchunk++;
    index = freehead[ibin];
    int ipage = index/chunkperpage;
    int ientry = index % chunkperpage;
    freehead[ibin] = freelist[index];
    return &pages[ipage][ientry*chunksize[ibin]];
  }

  // return indexed chunk to pool via free list
  // index = -1 if no allocated chunk

  void put(int index) {
    if (index < 0) return;
    int ipage = index/chunkperpage;
    int ibin = whichbin[ipage];
    nchunk--;
    ndatum -= chunksize[ibin];
    freelist[index] = freehead[ibin];
    freehead[ibin] = index;
  }

 private:
  int minchunk;       // min # of datums per chunk
  int maxchunk;       // max # of datums per chunk
  int nbin;           // # of bins to split min-to-max into
  int chunkperpage;   // # of chunks on every page, regardless of which bin
  int pagedelta;      // # of pages to allocate at once, default = 1
  int binsize;        // delta in chunk sizes between adjacent bins

  T **pages;          // list of allocated pages
  int *whichbin;      // which bin each page belongs to
  int npage;          // # of allocated pages
  int *freelist;      // each chunk points to next unused chunk in same bin
  int *freehead;      // index of first unused chunk in each bin
  int *chunksize;     // size of chunks in each bin

  void allocate(int ibin) {
    int oldpage = npage;
    npage += pagedelta;
    freelist = (int *) realloc(freelist,npage*chunkperpage*sizeof(int));
    pages = (T **) realloc(pages,npage*sizeof(T *));
    whichbin = (int *) realloc(whichbin,npage*sizeof(int));
    if (!freelist || !pages) {
      errorflag = 2;
      return;
    }

    // allocate pages with appropriate chunksize for ibin

    for (int i = oldpage; i < npage; i++) {
      whichbin[i] = ibin;
#if defined(LAMMPS_MEMALIGN)
      void *ptr;
      if (posix_memalign(&ptr, LAMMPS_MEMALIGN,
                         chunkperpage*chunksize[ibin]*sizeof(T)))
        errorflag = 2;
      pages[i] = (T *) ptr;
#else
      pages[i] = (T *) malloc(chunkperpage*chunksize[ibin]*sizeof(T));
      size += chunkperpage*chunksize[ibin];
      if (!pages[i]) errorflag = 2;
#endif
    }

    // reset free list for unused chunks on new pages

    freehead[ibin] = oldpage*chunkperpage;
    for (int i = freehead[ibin]; i < npage*chunkperpage; i++) freelist[i] = i+1;
    freelist[npage*chunkperpage-1] = -1;
  }
};

}

#endif
