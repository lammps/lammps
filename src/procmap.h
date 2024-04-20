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

#ifndef LMP_PROCMAP_H
#define LMP_PROCMAP_H

#include "pointers.h"

namespace LAMMPS_NS {

class ProcMap : protected Pointers {
 public:
  ProcMap(class LAMMPS *);

  void onelevel_grid(int, int *, int *, int, int, int *, int *);
  void twolevel_grid(int, int *, int *, int, int *, int *, int, int, int *, int *);
  void numa_grid(int, int, int *, int *, int *);
  void custom_grid(char *, int, int *, int *);
  void cart_map(int, int *, int *, int[3][2], int ***);
  void cart_map(int, int *, int, int *, int *, int[3][2], int ***);
  void xyz_map(char *, int *, int *, int[3][2], int ***);
  void xyz_map(char *, int *, int, int *, int *, int[3][2], int ***);
  void numa_map(int, int *, int *, int[3][2], int ***);
  void custom_map(int *, int *, int[3][2], int ***);
  void output(char *, int *, int ***);

 private:
  int procs_per_node;    // NUMA params
  int procs_per_numa;
  int node_id;        // which node I am in
  int nodegrid[3];    // 3d grid of nodes

  int **cmap;    // info in custom grid file

  int factor(int, int **);
  int combine_factors(int, int **, int, int **, int **);
  int cull_2d(int, int **, int);
  int cull_user(int, int **, int, int *);
  int cull_other(int, int **, int, int, int *, int *);
  int best_factors(int, int **, int *, int, int, int);
  void grid_shift(int, int, int &, int &);
};

}    // namespace LAMMPS_NS

#endif
