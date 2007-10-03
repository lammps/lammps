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

#ifndef MEMORY_H
#define MEMORY_H

#include "pointers.h"

namespace LAMMPS_NS {

class Memory : protected Pointers {
 public:
  Memory(class LAMMPS *);

  void *smalloc(int n, const char *);
  void sfree(void *);
  void *srealloc(void *, int n, const char *name);

  double *create_1d_double_array(int, int, const char *);
  void destroy_1d_double_array(double *, int);
  
  double **create_2d_double_array(int, int, const char *);
  void destroy_2d_double_array(double **);
  double **grow_2d_double_array(double **, int, int, const char *);

  int **create_2d_int_array(int, int, const char *);
  void destroy_2d_int_array(int **);
  int **grow_2d_int_array(int **, int, int, const char *);

  double **create_2d_double_array(int, int, int, const char *);
  void destroy_2d_double_array(double **, int);

  double ***create_3d_double_array(int, int, int, const char *);
  void destroy_3d_double_array(double ***);
  double ***grow_3d_double_array(double ***, int, int, int, const char *);

  double ***create_3d_double_array(int, int, int, int, const char *);
  void destroy_3d_double_array(double ***, int);

  double ***create_3d_double_array(int, int, int, int, int, int, const char *);
  void destroy_3d_double_array(double ***, int, int, int);

  int ***create_3d_int_array(int, int, int, const char *);
  void destroy_3d_int_array(int ***);

  double ****create_4d_double_array(int, int, int, int, const char *);
  void destroy_4d_double_array(double ****);
};

}

#endif
