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
  double cutoff_custom;        // cutoff set by requestor  

  // Variables for NBinStandard

  int nbinx,nbiny,nbinz;           // # of global bins
  int mbins;                       // # of local bins and offset on this proc
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  double binsizex,binsizey,binsizez;  // bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  int *binhead;                // index of first atom in each bin
  int *bins;                   // index of next atom in same bin
  int *atom2bin;               // bin assignment for each atom (local+ghost)

  // Analogues for NBinMultimulti2
  
  int *nbinx_multi2, *nbiny_multi2, *nbinz_multi2;
  int *mbins_multi2;
  int *mbinx_multi2, *mbiny_multi2, *mbinz_multi2;
  int *mbinxlo_multi2, *mbinylo_multi2, *mbinzlo_multi2;
  double *binsizex_multi2, *binsizey_multi2, *binsizez_multi2;
  double *bininvx_multi2, *bininvy_multi2, *bininvz_multi2;

  int **binhead_multi2;
  int **bins_multi2;
  int **atom2bin_multi2;
  
  NBin(class LAMMPS *);
  ~NBin();
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();

  virtual void bin_atoms_setup(int) = 0;
  virtual void setup_bins(int) = 0;
  virtual void bin_atoms() = 0;
  virtual double memory_usage() {return 0.0;}

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

  // data for standard NBin
  
  int maxatom;                      // size of bins array

  // data for standard NBin

  int maxbin;                       // size of binhead array

  // data for multi/multi2 NBin

  int maxtypes;                     // size of multi2 arrays
  int * maxbins_multi2;             // size of 2nd dimension of binhead_multi2 array

  // methods

  int coord2bin(double *);
  int coord2bin_multi2(double *, int);  
};

}

#endif

/* ERROR/WARNING messages:

E: Non-numeric positions - simulation unstable

UNDOCUMENTED

*/
