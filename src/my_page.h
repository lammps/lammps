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

/* ----------------------------------------------------------------------
MyPage = templated class for storing chunks of datums in pages
  chunks are not returnable, can only reset and start over
  replaces many small mallocs with a few large mallocs
  pages are never freed, so can reuse w/out reallocs
usage:
  request one datum at a time, repeat, clear
  request chunks of datums in each get() or vget(), repeat, clear
  chunk size can vary from request to request
  chunk size can be known in advance or registered after usage via vgot()
inputs:
   template T = one datum, e.g. int, double, struct, int[3]
     for int[3], access datum as ivec[i][2]
   maxchunk = max # of datums in one chunk, default = 1
   pagesize = # of datums in one page, default = 1024
              should be insure big enough to store multiple chunks
   pagedelta = # of pages to allocate at a time, default = 1
methods:
   T *get() = return ptr to one datum
   T *get(N) = return ptr to N datums, N < maxchunk required
   T *vget() = return ptr to maxchunk datums, use as needed, then call vgot()
     all gets return NULL if error encountered
   vgot(N) = used N datums of previous vget(), N < maxchunk required
   void reset() = clear pages w/out freeing
   int size() = return total size of allocated pages in bytes
public variables:
   ndatum = total # of stored datums
   nchunk = total # of stored chunks
   errorflag = flag for various error conditions
------------------------------------------------------------------------- */

#ifndef LAMMPS_MY_PAGE_H
#define LAMMPS_MY_PAGE_H

#include "stdlib.h"
namespace LAMMPS_NS {

template<class T>
class MyPage {
 public:
  int ndatum;      // total # of stored datums
  int nchunk;      // total # of stored chunks
  int errorflag;   // flag > 1 if error has occurred
                   // 1 = invalid inputs
                   // 2 = memory allocation error
                   // 3 = chunk size exceeded maxchunk

  MyPage(int user_maxchunk = 1, int user_pagesize = 1024, 
         int user_pagedelta = 1) {
    maxchunk = user_maxchunk;
    pagesize = user_pagesize;
    pagedelta = user_pagedelta;

    errorflag = 0;
    if (maxchunk <= 0 || pagesize <= 0 || pagedelta <= 0) errorflag = 1;
    if (maxchunk > pagesize) errorflag = 1;
    if (errorflag) return;

    // initial page allocation

    ndatum = nchunk = 0;
    pages = NULL;
    npage = 0;
    allocate();
    if (errorflag) return;
    ipage = index = 0;
    page = pages[ipage];
  }

  // free all pages of allocated memory

  ~MyPage() {
    for (int i = 0; i < npage; i++) free(pages[i]);
    free(pages);
  }

  // get ptr to one datum
  // return NULL if run out of memory

  T *get() {
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

  // get ptr to location that can store N datums
  // error if N > maxchunk
  // return NULL if run out of memory

  T *get(int n) {
    if (n > maxchunk) {
      errorflag = 3;
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

  // get ptr to location that can store maxchunk datums
  // will return same ptr as previous call if vgot() not called
  // return NULL if run out of memory
  
  T *vget() {
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

  // increment by N = # of values stored in loc returned by vget()
  // OK to not call if vget() ptr was not used
  // error if N > maxchunk

  void vgot(int n) {
    if (n > maxchunk) errorflag = 3;
    ndatum += n;
    nchunk++;
    index += n;
  }

  // reset index to beginning of first page
  // effectively clears all pages, without freeing any memory

  void reset() {
    ndatum = nchunk = 0;
    index = ipage = 0;
    page = pages[ipage];
  }

  // return total size of all allocated pages

  int size() {
    return npage*pagesize*sizeof(T);
  }

 private:
  int maxchunk;   // max # of datums in one requested chunk
  int pagesize;   // # of datums in one page, default = 1024
  int pagedelta;  // # of pages to allocate at once, default = 1

  T **pages;      // list of allocated pages
  T *page;        // ptr to current page
  int npage;      // # of allocated pages
  int ipage;      // index of current page
  int index;      // current index on current page
  
  void allocate() {
    npage += pagedelta;
    pages = (T **) realloc(pages,npage*sizeof(T *));
    if (!pages) {
      errorflag = 2;
      return;
    }
    for (int i = npage-pagedelta; i < npage; i++) {
      pages[i] = (T *) malloc(pagesize*sizeof(T));
      if (!pages[i]) errorflag = 2;
    }
  }
};

}

#endif
