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

#ifndef LMP_NEIGH_BIN_H
#define LMP_NEIGH_BIN_H

#include "pointers.h"

namespace LAMMPS_NS {

class NeighBin : protected Pointers {
 public:
  int style;                       // ID for NeighBin method this is

  bigint last_setup,last_bin;      // timesteps for last operations performed
  bigint last_bin_memory;

  int nbinx,nbiny,nbinz;           // # of global bins
  int mbins;                       // # of local bins and offset on this proc
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  double binsizex,binsizey,binsizez;  // bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  int *binhead;                    // index of first atom in each bin
  int *bins;                       // index of next atom in same bin

  NeighBin(class LAMMPS *);
  ~NeighBin();
  void copy_neighbor_info();
  virtual void bin_atoms_setup(int);
  virtual bigint memory_usage();

  virtual void setup_bins(int);
  virtual void bin_atoms() = 0;

 protected:

  // data from Neighbor class

  int includegroup;
  double cutneighmin;
  double cutneighmax;
  int binsizeflag;
  double binsize_user;
  double *bboxlo,*bboxhi;

  // data common to all NeighBin variants

  int dimension;
  int triclinic;

  int maxbin;                       // size of binhead array
  int maxatom;                      // size of bins array

  // methods

  int coord2bin(double *);
};

}

#endif

/* ERROR/WARNING messages:

*/
