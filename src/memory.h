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

#include "pointers.h"

namespace LAMMPS_NS {

class Memory : protected Pointers {
 public:
  Memory(class LAMMPS *);

  // older routines

  void *smalloc(int n, const char *);
  void sfree(void *);
  void *srealloc(void *, int n, const char *);

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

  // newer routines

 public:
  void *smalloc_new(bigint n, const char *);
  void *srealloc_new(void *, bigint n, const char *);
  void sfree_new(void *);


/* ----------------------------------------------------------------------
   to avoid code bloat, only use these for int,double,float,char
   not for int* or double** or arbitrary structs
   in those cases, just use smalloc,srealloc,sfree directly
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   create a 1d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name)
    {
      bigint nbytes = sizeof(TYPE) * n;
      array = (TYPE *) smalloc_new(nbytes,name);
      return array;
    };

/* ----------------------------------------------------------------------
   grow or shrink 1d array
   if dim is 0, return NULL 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *grow(TYPE *&array, int n, const char *name) 
    {
      if (array == NULL) return create(array,n,name);
      if (n == 0) {
	destroy(array);
	return NULL;
      }
      
      bigint nbytes = sizeof(TYPE) * n;
      array = (TYPE *) srealloc_new(array,nbytes,name);
      return array;
    };

/* ----------------------------------------------------------------------
   destroy a 1d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE *array) 
    {
      sfree_new(array);
    };

/* ----------------------------------------------------------------------
   create a 1d array with index from nlo to nhi inclusive 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *create(TYPE *&array, int nlo, int nhi, const char *name) 
    {
      bigint nbytes = sizeof(TYPE) * (nhi-nlo+1);
      array = (TYPE *) smalloc_new(nbytes,name);
      array = array-nlo;
      return array;
    }

/* ----------------------------------------------------------------------
   destroy a 1d array with index offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE *array, int offset) 
    {
      if (array) sfree_new(array+offset);
    }

/* ----------------------------------------------------------------------
   create a 2d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name) 
    {
      bigint nbytes = sizeof(TYPE) * n1*n2;
      TYPE *data = (TYPE *) smalloc_new(nbytes,name);
      nbytes = sizeof(TYPE *) * n1;
      array = (TYPE **) smalloc_new(nbytes,name);
      
      int n = 0;
      for (int i = 0; i < n1; i++) {
	array[i] = &data[n];
	n += n2;
      }
      return array;
    }

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
   if either dim is 0, destroy it and return NULL 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **grow(TYPE **&array, int n1, int n2, const char *name) 
    {
      if (array == NULL) return create(array,n1,n2,name);
      if (n1 == 0 || n2 == 0) {
	destroy(array);
	return NULL;
      }
      
      bigint nbytes = sizeof(TYPE) * n1*n2;
      TYPE *data = (TYPE *) srealloc_new(array[0],nbytes,name);
      nbytes = sizeof(TYPE *) * n1;
      array = (TYPE **) srealloc_new(array,nbytes,name);
      
      int n = 0;
      for (int i = 0; i < n1; i++) {
	array[i] = &data[n];
	n += n2;
      }
      return array;
    }

/* ----------------------------------------------------------------------
   destroy a 2d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE **array)
    {
      if (array == NULL) return;
      sfree_new(array[0]);
      sfree_new(array);
    }

/* ----------------------------------------------------------------------
   create a 2d array with 2nd index from n2lo to n2hi inclusive 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2lo, int n2hi,
		  const char *name)
    {
      int n2 = n2hi - n2lo + 1;
      create(array,n1,n2,name);
      for (int i = 0; i < n1; i++) array[i] -= n2lo;
      return array;
    }

/* ----------------------------------------------------------------------
   destroy a 2d array with 2nd index offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE **array, int offset)
    {
      if (array == NULL) return;
      sfree_new(&array[0][offset]);
      sfree_new(array);
    }

/* ----------------------------------------------------------------------
   create a 3d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, int n3, const char *name) 
    {
      bigint nbytes = sizeof(TYPE) * n1*n2*n3;
      TYPE *data = (TYPE *) smalloc_new(nbytes,name);
      nbytes = sizeof(TYPE *) * n1*n2;
      TYPE **plane = (TYPE **) smalloc_new(nbytes,name);
      nbytes = sizeof(TYPE **) * n1;
      array = (TYPE ***) smalloc_new(nbytes,name);
      
      int i,j;
      int n = 0;
      for (i = 0; i < n1; i++) {
	array[i] = &plane[i*n2];
	for (j = 0; j < n2; j++) {
	  plane[i*n2+j] = &data[n];
	  n += n3;
	}
      }
      return array;
    }

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 3d array
   last 2 dims must stay the same
   if any dim is 0, destroy it and return NULL 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***grow(TYPE ***&array, int n1, int n2, int n3, const char *name) 
    {
      if (array == NULL) return create(array,n1,n2,n3,name);
      if (n1 == 0 || n2 == 0 || n3 == 0) {
	destroy(array);
	return NULL;
      }
      
      bigint nbytes = sizeof(TYPE) * n1*n2*n3;
      TYPE *data = (TYPE *) srealloc_new(array[0][0],nbytes,name);
      nbytes = sizeof(TYPE *) * n1*n2;
      TYPE **plane = (TYPE **) srealloc_new(array[0],nbytes,name);
      nbytes = sizeof(TYPE **) * n1;
      array = (TYPE ***) srealloc_new(array,nbytes,name);
      
      int i,j;
      int n = 0;
      for (i = 0; i < n1; i++) {
	array[i] = &plane[i*n2];
	for (j = 0; j < n2; j++) {
	  plane[i*n2+j] = &data[n];
	  n += n3;
	}
      }
      return array;
    }

/* ----------------------------------------------------------------------
   destroy a 3d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ***array) 
    {
      if (array == NULL) return;
      sfree_new(array[0][0]);
      sfree_new(array[0]);
      sfree_new(array);
    }

/* ----------------------------------------------------------------------
   create a 3d array with 1st index from n1lo to n1hi inclusive 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1lo, int n1hi, 
		   int n2, int n3, const char *name)
    {
      int n1 = n1hi - n1lo + 1;
      create(array,n1,n2,n3,name);
      return array-n1lo;
    }

/* ----------------------------------------------------------------------
   free a 3d array with 1st index offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ***array, int offset)
    {
      if (array) destroy(array+offset);
    }

/* ----------------------------------------------------------------------
   create a 3d array with
   1st index from n1lo to n1hi inclusive,
   2nd index from n2lo to n2hi inclusive,
   3rd index from n3lo to n3hi inclusive 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1lo, int n1hi, int n2lo, int n2hi,
		   int n3lo, int n3hi, const char *name)
    {
      int n1 = n1hi - n1lo + 1;
      int n2 = n2hi - n2lo + 1;
      int n3 = n3hi - n3lo + 1;
      create(array,n1,n2,n3,name);
      
      for (int i = 0; i < n1*n2; i++) array[0][i] -= n3lo;
      for (int i = 0; i < n1; i++) array[i] -= n2lo;
      return array-n1lo;
    }

/* ----------------------------------------------------------------------
   free a 3d array with all 3 indices offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ***array, int n1_offset, int n2_offset, int n3_offset)
    {
      if (array == NULL) return;
      sfree_new(&array[n1_offset][n2_offset][n3_offset]);
      sfree_new(&array[n1_offset][n2_offset]);
      sfree_new(array + n1_offset);
    }

/* ----------------------------------------------------------------------
   create a 4d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ****create(TYPE ****&array, int n1, int n2, int n3, int n4,
		    const char *name)
    {
      bigint nbytes = sizeof(TYPE) * n1*n2*n3*n4;
      TYPE *data = (double *) smalloc_new(nbytes,name);
      nbytes = sizeof(TYPE *) * n1*n2*n3;
      TYPE **cube = (double **) smalloc_new(nbytes,name);
      nbytes = sizeof(TYPE **) * n1*n2;
      TYPE ***plane = (double ***) smalloc_new(nbytes,name);
      nbytes = sizeof(TYPE ***) * n1;
      array = (double ****) smalloc_new(nbytes,name);
      
      int i,j,k;
      int n = 0;
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
      return array;
    }

/* ----------------------------------------------------------------------
   destroy a 4d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ****array)
    {
      if (array == NULL) return;
      sfree_new(array[0][0][0]);
      sfree_new(array[0][0]);
      sfree_new(array[0]);
      sfree_new(array);
    }

/* ----------------------------------------------------------------------
   memory usage of arrays, including pointers
------------------------------------------------------------------------- */

  template <typename TYPE>
    bigint usage(TYPE *array, int n)
    {
      bigint bytes = sizeof(TYPE) * n;
      return bytes;
    }

  template <typename TYPE>
    bigint usage(TYPE **array, int n1, int n2)
    {
      bigint bytes = sizeof(TYPE) * n1*n2;
      bytes += sizeof(TYPE *) * n1;
      return bytes;
    }

  template <typename TYPE>
    bigint usage(TYPE ***array, int n1, int n2, int n3)
    {
      bigint bytes = sizeof(TYPE) * n1*n2*n3;
      bytes += sizeof(TYPE *) * n1*n2;
      bytes += sizeof(TYPE **) * n1;
      return bytes;
    }

  template <typename TYPE>
    bigint usage(TYPE ****array, int n1, int n2, int n3, int n4)
    {
      bigint bytes = sizeof(TYPE) * n1*n2*n3*n4;
      bytes += sizeof(TYPE *) * n1*n2*n3;
      bytes += sizeof(TYPE **) * n1*n2;
      bytes += sizeof(TYPE ***) * n1;
      return bytes;
    }
};

}

#endif
