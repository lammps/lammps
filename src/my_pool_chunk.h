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

#ifndef LAMMPS_MY_POOL_CHUNK_H
#define LAMMPS_MY_POOL_CHUNK_H

namespace LAMMPS_NS {

template <class T> class MyPoolChunk {
 public:
  int ndatum;    // total # of stored datums
  int nchunk;    // total # of stored chunks

  MyPoolChunk(int user_minchunk = 1, int user_maxchunk = 1, int user_nbin = 1,
              int user_chunkperpage = 1024, int user_pagedelta = 1);

  // free all allocated memory

  ~MyPoolChunk();

  // return pointer/index of unused chunk of size maxchunk

  T *get(int &index);

  // return pointer/index of unused chunk of size N

  T *get(int n, int &index);

  // return indexed chunk to pool via free list
  // index = -1 if no allocated chunk

  void put(int index);

  // total memory used in bytes

  double size() const;

  /** Return error status
   *
   * \return 0 if no error, 1 if invalid input, 2 if malloc() failed, 3 if chunk > maxchunk */

  int status() const { return errorflag; }

 private:
  int minchunk;        // min # of datums per chunk
  int maxchunk;        // max # of datums per chunk
  int nbin;            // # of bins to split min-to-max into
  int chunkperpage;    // # of chunks on every page, regardless of which bin
  int pagedelta;       // # of pages to allocate at once, default = 1
  int binsize;         // delta in chunk sizes between adjacent bins
  int errorflag;       // flag > 0 if error has occurred
                       // 1 = invalid inputs
                       // 2 = memory allocation error
                       // 3 = chunk size exceeded maxchunk

  T **pages;         // list of allocated pages
  int *whichbin;     // which bin each page belongs to
  int npage;         // # of allocated pages
  int *freelist;     // each chunk points to next unused chunk in same bin
  int *freehead;     // index of first unused chunk in each bin
  int *chunksize;    // size of chunks in each bin

  void allocate(int ibin);
};
}    // namespace LAMMPS_NS
#endif
