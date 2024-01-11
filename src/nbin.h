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

#ifndef LMP_NBIN_H
#define LMP_NBIN_H

#include "pointers.h"    // IWYU pragma: keep

namespace LAMMPS_NS {

class NBin : protected Pointers {
 public:
  int istyle;              // 1-N index into binnames
  bigint last_bin;         // last timestep atoms were binned
  double cutoff_custom;    // cutoff set by requestor

  // Variables for NBinStandard

  int nbinx, nbiny, nbinz;    // # of global bins
  int mbins;                  // # of local bins and offset on this proc
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;

  double binsizex, binsizey, binsizez;    // bin sizes and inverse sizes
  double bininvx, bininvy, bininvz;

  int *binhead;     // index of first atom in each bin
  int *bins;        // index of next atom in same bin
  int *atom2bin;    // bin assignment for each atom (local+ghost)

  // Analogues for NBinMultimulti

  int *nbinx_multi, *nbiny_multi, *nbinz_multi;
  int *mbins_multi;
  int *mbinx_multi, *mbiny_multi, *mbinz_multi;
  int *mbinxlo_multi, *mbinylo_multi, *mbinzlo_multi;
  double *binsizex_multi, *binsizey_multi, *binsizez_multi;
  double *bininvx_multi, *bininvy_multi, *bininvz_multi;

  int **binhead_multi;

  NBin(class LAMMPS *);
  ~NBin() override;
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();

  virtual void bin_atoms_setup(int) = 0;
  virtual void setup_bins(int) = 0;
  virtual void bin_atoms() = 0;
  virtual double memory_usage() { return 0.0; }

  // Kokkos package

  int kokkos;    // 1 if class stores Kokkos data

 protected:
  // data from Neighbor class

  int includegroup;
  double cutneighmin;
  double cutneighmax;
  int binsizeflag;
  double binsize_user;
  double *bboxlo, *bboxhi;
  int ncollections;
  double **cutcollectionsq;

  // data common to all NBin variants

  int dimension;
  int triclinic;

  // data for standard NBin

  int maxatom;    // size of bins array
  int maxbin;     // size of binhead array

  // data for multi NBin

  int maxcollections;    // size of multi arrays
  int *maxbins_multi;    // size of 2nd dimension of binhead_multi array

  // methods

  int coord2bin(double *);
  int coord2bin_multi(double *, int);
};

}    // namespace LAMMPS_NS

#endif
