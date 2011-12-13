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

#ifndef LMP_PROCMAP_H
#define LMP_PROCMAP_H

#include "pointers.h"

namespace LAMMPS_NS {

class ProcMap : protected Pointers {
 public:
  ProcMap(class LAMMPS *);
  ~ProcMap() {}
  int onelevel_grid(int, int *, int *, int, int, int *);
  int twolevel_grid(int, int *, int *, int, int *, int *, int, int, int *);
  int numa_grid(int, int *, int *, int *);
  void custom_grid(char *, int, int *, int *);
  void cart_map(int, int *, int *, int [3][2], int ***);
  void cart_map(int, int *, int *, int *, int [3][2], int ***);
  void xyz_map(char *, int *, int *, int [3][2], int ***);
  void xyz_map(char *, int *, int *, int *, int [3][2], int ***);
  void numa_map(int *, int *, int [3][2], int ***);
  void custom_map(int *, int [3][2], int ***);
  void output(int ***, char *);

 private:
  int other_style;
  int other_procgrid[3];
  int nodegrid[3];
  int node_id;
  int procs_per_node;
  int procs_per_numa;

  int procs2box(int, int *, int *, const int, const int, const int, int);
  void grid_shift(int, int, int &, int &);
};

}

#endif
