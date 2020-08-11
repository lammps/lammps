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

#ifndef LMP_NBIN_H
#define LMP_NBIN_H

#include "pointers.h"

namespace LAMMPS_NS {

class NBin : protected Pointers {
 public:
  int istyle;                      // 1-N index into binnames
  bigint last_bin;                 // last timestep atoms were binned

  int nbinx,nbiny,nbinz;           // # of global bins
  int mbins;                       // # of local bins and offset on this proc
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  double binsizex,binsizey,binsizez;  // bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  int *binhead;                // index of first atom in each bin
  int *bins;                   // index of next atom in same bin
  int *atom2bin;               // bin assignment for each atom (local+ghost)

  //CAC package bin arrays
  int *bin_ncontent;          //number of contents in each bin
  int **bin_content;          //array of local and ghost indices in each bin
  int *quad2bin;              //bin location of each local quadrature point
  int *nbin_element_overlap;  //array storing the number of bins this element overlaps
  int **bin_element_overlap;  //set of bins this element overlaps

  double cutoff_custom;        // cutoff set by requestor

  NBin(class LAMMPS *);
  virtual ~NBin();
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();
  virtual void bin_atoms_setup(int);
  bigint memory_usage();

  virtual void setup_bins(int) = 0;
  virtual void bin_atoms() = 0;
  virtual int coord2bin(double *);

  // Kokkos package

  int kokkos;                       // 1 if class stores Kokkos data

 protected:

  // data from Neighbor class

  int includegroup;
  double cutneighmin;
  double cutneighmax;
  int binsizeflag;
  double binsize_user;
  double *bboxlo,*bboxhi;

  // data common to all NBin variants

  int dimension;
  int triclinic;

  int maxbin;                       // size of binhead array
  int maxatom;                      // size of bins array
};

}

#endif

/* ERROR/WARNING messages:

E: Non-numeric positions - simulation unstable

UNDOCUMENTED

*/
