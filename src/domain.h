/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef DOMAIN_H
#define DOMAIN_H

#include "lammps.h"
class Lattice;
class Region;

class Domain : public LAMMPS {
 public:
  int box_exist;                            // 0 = not yet created, 1 = exists

  int nonperiodic;                          // 0 = periodic in all 3 dims
                                            // 1 = periodic or fixed in all 6
                                            // 2 = shrink-wrap in any of 6
  int xperiodic,yperiodic,zperiodic;        // 0 = not periodic, 1 = periodic
  int boundary[3][2];                       // settings for 6 boundaries
                                            // 0 = periodic
                                            // 1 = fixed non-periodic
                                            // 2 = shrink-wrap non-periodic
                                            // 3 = shrink-wrap non-per w/ min

  double minxlo,minxhi;                     // minimum size of global box
  double minylo,minyhi;                     // when shrink-wrapping
  double minzlo,minzhi;

  double boxxlo,boxxhi;                     // global box boundaries 
  double boxylo,boxyhi;
  double boxzlo,boxzhi;

  double xprd,yprd,zprd;                    // global box size
  double xprd_half,yprd_half,zprd_half;

  double subxlo,subxhi;                     // sub-box boudaries on this proc
  double subylo,subyhi;
  double subzlo,subzhi;

  double boxlo[3],boxhi[3];                 // global box bounds as arrays
  double prd[3];                            // global box size as array
  double sublo[3],subhi[3];                 // sub-box bounds as arrays
  int periodicity[3];                       // xyz periodic as array

  int box_change;            // 1 if box bounds ever change, 0 if fixed

  Lattice *lattice;              // user-defined lattice

  int nregion;                   // # of defined Regions
  int maxregion;                 // max # list can hold
  Region **regions;              // list of defined Regions

  Domain();
  ~Domain();
  void init();
  void set_initial_box();
  void set_global_box();
  void set_local_box();
  void reset_box();
  void pbc();
  void remap(double &, double &, double &, int &);
  void unmap(double &, double &, double &, int);
  void minimum_image(double *, double *, double *);
  void minimum_image(double *);
  void set_lattice(int, char **);
  void add_region(int, char **);
  void set_boundary(int, char **);
};

#endif
