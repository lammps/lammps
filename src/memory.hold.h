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

#ifndef MEMORY_H
#define MEMORY_H

#include "lammps.h"

class Memory : public LAMMPS {
 public:
  Memory() {}
  ~Memory() {}

  void *smalloc(int n, char *);
  void sfree(void *);
  void *srealloc(void *, int n, char *name);

  double *create_1d_double_array(int, int, char *);
  void destroy_1d_double_array(double *, int);
  
  double **create_2d_double_array(int, int, char *);
  void destroy_2d_double_array(double **);
  double **grow_2d_double_array(double **, int, int, char *);

  int **create_2d_int_array(int, int, char *);
  void destroy_2d_int_array(int **);
  int **grow_2d_int_array(int **, int, int, char *);

  double **create_2d_double_array(int, int, int, char *);
  void destroy_2d_double_array(double **, int);

  double ***create_3d_double_array(int, int, int, char *);
  void destroy_3d_double_array(double ***);
  double ***grow_3d_double_array(double ***, int, int, int, char *);

  double ***create_3d_double_array(int, int, int, int, char *);
  void destroy_3d_double_array(double ***, int);

  double ***create_3d_double_array(int, int, int, int, int, int, char *);
  void destroy_3d_double_array(double ***, int, int, int);
};

#endif
