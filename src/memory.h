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

#ifndef LMP_MEMORY_H
#define LMP_MEMORY_H

#include "lmptype.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Memory : protected Pointers {
 public:
  Memory(class LAMMPS *);

  void *smalloc(bigint n, const char *);
  void sfree(void *);
  void *srealloc(void *, bigint n, const char *);

  /* templated allocators and deallocators */

  /// create a 1d array with index from nlo to nhi inclusive 
  template <typename T>
  bool create_1d_array(T **ptr, bigint nlo, bigint nhi, const char *name) {
    if (ptr == 0) return false;

    bigint n = nhi - nlo + 1;
    T *array = (T *) smalloc(n*sizeof(T),name);
    *ptr=array-nlo;
    return true;
  };

  /// free 1d array with index offset.
  template <typename T>
  void destroy_1d_array(T *array, bigint offset) {
    if (array == NULL) return;
    sfree(array+offset);
  };

  /// create a 2d array. if either dim is 0, fail.
  template <typename T>
  bool create_2d_array(T ***ptr, bigint n1, bigint n2, const char *name) {

    if (ptr == 0) return false;

    if (n1 == 0 || n2 == 0) {
      *ptr = NULL;
      return false;
    }

    T *data = (T *) smalloc(n1*n2*sizeof(T),name);
    T **array = (T **) smalloc(n1*sizeof(T *),name);

    bigint n = 0;
    for (bigint i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    *ptr=array;
    return true;
  };

  /// grow or shrink 1st dim of a 2d array. last dim must stay the same.
  template <typename T>
  bool grow_2d_array(T ***array, bigint n1, bigint n2, const char *name) {

    if (array == 0) return false;

    if (n1 == 0 || n2 == 0) {
      destroy_2d_array(*array);
      *array = NULL;
      return false;
    }

    if (*array == NULL) return create_2d_array(array,n1,n2,name);

    T *data = (T *) srealloc((*array)[0],n1*n2*sizeof(T),name);
    *array = (T **) srealloc(*array,n1*sizeof(T *),name);

    bigint n = 0;
    for (bigint i = 0; i < n1; i++) {
      (*array)[i] = &data[n];
      n += n2;
    }
    return true;
  };

  /// free regular 2d array.
  template <typename T>
  void destroy_2d_array(T **array) {
    if (array == NULL) return;
    sfree(array[0]);
    sfree(array);
  };

  /// create a 2d array with 2nd index from n2lo to n2hi inclusive 
  template <typename T>
  bool create_2d_array(T ***ptr, bigint n1, bigint n2lo, bigint n2hi,
		       const char *name) {

    if (ptr == 0) return false;

    bigint n2 = n2hi - n2lo + 1;
    T **array;
    if (create_2d_array(&array,n1,n2,name)) {
      for (bigint i = 0; i < n1; i++)
	array[i] -= n2lo;
      *ptr = array;
      return true;
    } else return false;
  };

  /// free a 2d array with 2nd index offset
  template <typename T>
  void destroy_2d_array(T **array, bigint offset) {
    if (array == NULL) return;
    sfree(&array[0][offset]);
    sfree(array);
  };

  /// create a 3d array 
  template <typename T>
  bool create_3d_array(T ****ptr, bigint n1, bigint n2, bigint n3,
		       const char *name) {

    if (ptr == 0) return false;

    if (n1 == 0 || n2 == 0 || n3 == 0) {
      *ptr = NULL;
      return false;
    }

    T *data = (T *) smalloc(n1*n2*n3*sizeof(T),name);
    T **plane = (T **) smalloc(n1*n2*sizeof(T *),name);
    T ***array = (T ***) smalloc(n1*sizeof(T **),name);

    bigint i,j;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      array[i] = &plane[i*n2];
      for (j = 0; j < n2; j++) {
	plane[i*n2+j] = &data[n];
	n += n3;
      }
    }
    *ptr = array;
    return true;
  };

  /// free a 3d double array 
  template <typename T>
  void destroy_3d_array(T ***array) {
    if (array == NULL) return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
  };

  /// grow or shrink 1st dim of a 3d array, last 2 dims must stay the same
  template <typename T>
  bool grow_3d_array(T ****array, bigint n1, bigint n2,
		     bigint n3, const char *name) {

    if (array == 0) return false;

    if (n1 == 0 || n2 == 0 || n3 == 0) {
      destroy_3d_array(*array);
      *array = NULL;
      return false;
    }

    if (*array == NULL) return create_3d_array(array,n1,n2,n3,name);

    T *data = (T *) srealloc((*array)[0][0],n1*n2*n3*sizeof(T),name);
    T **plane = (T **) srealloc((*array)[0],n1*n2*sizeof(T *),name);
    *array = (T ***) srealloc(*array,n1*sizeof(T **),name);

    bigint i,j;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      (*array)[i] = &plane[i*n2];
      for (j = 0; j < n2; j++) {
	plane[i*n2+j] = &data[n];
	n += n3;
      }
    }
    return true;
  };

  /// a 3d array with 1st index from n1lo to n1hi inclusive 
  template <typename T>
  bool create_3d_array(T ****ptr, bigint n1lo, bigint n1hi,
		       bigint n2, bigint n3, const char *name) {

    if (ptr == 0) return false;

    bigint n1 = n1hi - n1lo + 1;
    T ***array;
    if (create_3d_array(&array,n1,n2,n3,name)) {
      *ptr = array-n1lo;
      return true;
    } else return false;
  };

  /// free a 3d array with 1st index offset 
  template <typename T>
  void destroy_3d_array(T ***array, bigint offset) {
    if (array) destroy_3d_array(array + offset);
  };

  //! create a 3d array with 1st index from n1lo to n1hi inclusive,
  // 2nd index from n2lo to n2hi inclusive, 
  // 3rd index from n3lo to n3hi inclusive 
  template <typename T>
  bool create_3d_array(T ****ptr, bigint n1lo, bigint n1hi,
		       bigint n2lo, bigint n2hi, bigint n3lo,
		       bigint n3hi, const char *name) {

    if (ptr == 0) return false;

    bigint n1 = n1hi - n1lo + 1;
    bigint n2 = n2hi - n2lo + 1;
    bigint n3 = n3hi - n3lo + 1;
    T ***array;
    if (create_3d_array(&array,n1,n2,n3,name)) {
      bigint i;
      for (i = 0; i < n1*n2; i++) array[0][i] -= n3lo;
      for (i = 0; i < n1; i++) array[i] -= n2lo;
      *ptr = array-n1lo;
      return true;
    } else return false;
  };

  /// free a 3d array with all 3 indices offset 
  template <typename T>
  void destroy_3d_array(T ***array, bigint n1_offset,
			bigint n2_offset, bigint n3_offset) {
    if (array == NULL) return;
    sfree(&array[n1_offset][n2_offset][n3_offset]);
    sfree(&array[n1_offset][n2_offset]);
    sfree(array + n1_offset);
  };

  /// create a 4d array 
  template <typename T>
  bool create_4d_array(T *****ptr, bigint n1, bigint n2, bigint n3,
		       bigint n4, const char *name) {

    if (ptr == 0) return false;

    if (n1==0 || n2==0 || n3==0 || n4==0) {
      *ptr = NULL;
      return false;
    }

    T *data = (T *) smalloc(n1*n2*n3*n4*sizeof(T),name);
    T **cube = (T **) smalloc(n1*n2*n3*sizeof(T *),name);
    T ***plane = (T ***) smalloc(n1*n2*sizeof(T **),name);
    T ****array = (T ****) smalloc(n1*sizeof(T ***),name);

    bigint i,j,k;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      array[i] = &plane[i*n2];
      for (j = 0; j < n2; j++) {
	plane[i*n2+j] = &cube[i*n2*n3+j*n3];
	for (k = 0; k < n3; k++) {
	  cube[i*n2*n3+j*n3+k] = &data[n];
	  n += n4;
	}
      }
    }
    *ptr = array;
    return true;
  };

  /// free a 4d array 
  template <typename T>
  void destroy_4d_array(T ****array) {
    if (array == NULL) return;
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
  };

  /* regular allocators */
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

  float **create_2d_float_array(int, int, const char *);
  void destroy_2d_float_array(float **);
  float **grow_2d_float_array(float **, int, int, const char *);

  float **create_2d_float_array(int, int, int, const char *);
  void destroy_2d_float_array(float **, int);

  double ***create_3d_double_array(int, int, int, const char *);
  void destroy_3d_double_array(double ***);
  double ***grow_3d_double_array(double ***, int, int, int, const char *);

  double ***create_3d_double_array(int, int, int, int, const char *);
  void destroy_3d_double_array(double ***, int);

  double ***create_3d_double_array(int, int, int, int, int, int, const char *);
  void destroy_3d_double_array(double ***, int, int, int);

  float ***create_3d_float_array(int, int, int, const char *);
  void destroy_3d_float_array(float ***);
  float ***grow_3d_float_array(float ***, int, int, int, const char *);

  float ***create_3d_float_array(int, int, int, int, const char *);
  void destroy_3d_float_array(float ***, int);

  float ***create_3d_float_array(int, int, int, int, int, int, const char *);
  void destroy_3d_float_array(float ***, int, int, int);

  int ***create_3d_int_array(int, int, int, const char *);
  void destroy_3d_int_array(int ***);

  double ****create_4d_double_array(int, int, int, int, const char *);
  void destroy_4d_double_array(double ****);
};

}

#endif
