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

/* ----------------------------------------------------------------------
   templated class for storing chunks of datums in pages
------------------------------------------------------------------------- */

#ifndef LAMMPS_MY_PAGE_H
#define LAMMPS_MY_PAGE_H

#include "lmptype.h"

namespace LAMMPS_NS {

struct HyperOneCoeff {
  double biascoeff;
  tagint tag;
};

template <class T> class MyPage {
 public:
  int ndatum;    // total # of stored datums
  int nchunk;    // total # of stored chunks
  MyPage();
  virtual ~MyPage();

  int init(int user_maxchunk = 1, int user_pagesize = 1024, int user_pagedelta = 1);

  T *get(int n = 1);

  /** Get pointer to location that can store *maxchunk* items.
   *
   * This will return the same pointer as the previous call to
   * this function unless vgot() is called afterwards to record
   * how many items of the chunk were actually used.
   *
   * \return pointer to chunk of memory or null pointer if run out of memory */

  T *vget()
  {
    if (index + maxchunk <= pagesize) return &page[index];
    ipage++;
    if (ipage == npage) {
      allocate();
      if (errorflag) return nullptr;
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
   * \param  n  Number of items used in previously reserved chunk */

  void vgot(int n)
  {
    if (n > maxchunk) errorflag = 1;
    ndatum += n;
    nchunk++;
    index += n;
  }

  void reset();

  /** Return total size of allocated pages
   *
   * \return total storage used in bytes */

  double size() const { return (double) npage * pagesize * sizeof(T); }

  /** Return error status
   *
   * \return 0 if no error, 1 requested chunk size > maxchunk, 2 if malloc failed */

  int status() const { return errorflag; }

 private:
  T **pages;    // list of allocated pages
  T *page;      // ptr to current page
  int npage;    // # of allocated pages
  int ipage;    // index of current page
  int index;    // current index on current page

  int maxchunk;     // max # of datums in one requested chunk
  int pagesize;     // # of datums in one page, default = 1024
  int pagedelta;    // # of pages to allocate at once, default = 1

  int errorflag;    // flag > 0 if error has occurred
                    // 1 = chunk size exceeded maxchunk
                    // 2 = memory allocation error
  void allocate();
  void deallocate();
};

}    // namespace LAMMPS_NS

#endif
