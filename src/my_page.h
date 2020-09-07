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
   templated class for storing chunks of datums in pages
------------------------------------------------------------------------- */

#ifndef LAMMPS_MY_PAGE_H
#define LAMMPS_MY_PAGE_H

#if defined(LMP_USER_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif

namespace LAMMPS_NS {

template<class T>
class MyPage {
 public:
  int ndatum;      // total # of stored datums
  int nchunk;      // total # of stored chunks
  MyPage();
  virtual ~MyPage();

  int init(int user_maxchunk = 1, int user_pagesize = 1024,
           int user_pagedelta = 1);

  T *get();
  T *get(int n);

  T *vget();
  void vgot(int n);

  void reset();

  /** Return total size of allocated pages
   *
   * \return total storage used in bytes */

  int size() const {
    return npage*pagesize*sizeof(T);
  }

  /** Return error status
   *
   * \return 0 if no error, 1 requested chunk size > maxchunk, 2 if malloc failed */

  int status() const { return errorflag; }

 private:
  T **pages;      // list of allocated pages
  T *page;        // ptr to current page
  int npage;      // # of allocated pages
  int ipage;      // index of current page
  int index;      // current index on current page

  int maxchunk;   // max # of datums in one requested chunk
  int pagesize;   // # of datums in one page, default = 1024
  int pagedelta;  // # of pages to allocate at once, default = 1

  int errorflag;  // flag > 0 if error has occurred
                  // 1 = chunk size exceeded maxchunk
                  // 2 = memory allocation error
  void allocate();
};

}

#endif
