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

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Memory::Memory(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   safe malloc 
------------------------------------------------------------------------- */

void *Memory::smalloc(int n, const char *name)
{
  if (n == 0) return NULL;
  void *ptr = malloc(n);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate %d bytes for array %s",n,name);
    error->one(str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe free 
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   safe realloc 
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, int n, const char *name)
{
  if (n == 0) {
    sfree(ptr);
    return NULL;
  }

  ptr = realloc(ptr,n);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to reallocate %d bytes for array %s",n,name);
    error->one(str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   create a 1d double array with index from nlo to nhi inclusive 
------------------------------------------------------------------------- */

double *Memory::create_1d_double_array(int nlo, int nhi, const char *name)
{
  int n = nhi - nlo + 1;
  double *array = (double *) smalloc(n*sizeof(double),name);
  return array-nlo;
}

/* ----------------------------------------------------------------------
   free a 1d double array with index offset 
------------------------------------------------------------------------- */

void Memory::destroy_1d_double_array(double *array, int offset)
{
  if (array == NULL) return;
  sfree(array + offset);
}

/* ----------------------------------------------------------------------
   create a 2d double array 
------------------------------------------------------------------------- */

double **Memory::create_2d_double_array(int n1, int n2, const char *name)

{
  double *data = (double *) smalloc(n1*n2*sizeof(double),name);
  double **array = (double **) smalloc(n1*sizeof(double *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 2d double array 
------------------------------------------------------------------------- */

void Memory::destroy_2d_double_array(double **array)

{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d double array
   last dim must stay the same
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

double **Memory::grow_2d_double_array(double **array,
				      int n1, int n2, const char *name)

{
  if (array == NULL) return create_2d_double_array(n1,n2,name);

  double *data = (double *) srealloc(array[0],n1*n2*sizeof(double),name);
  array = (double **) srealloc(array,n1*sizeof(double *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   create a 2d int array
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

int **Memory::create_2d_int_array(int n1, int n2, const char *name)

{
  if (n1 == 0 || n2 == 0) return NULL;

  int *data = (int *) smalloc(n1*n2*sizeof(int),name);
  int **array = (int **) smalloc(n1*sizeof(int *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 2d int array 
------------------------------------------------------------------------- */

void Memory::destroy_2d_int_array(int **array)

{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d int array
   last dim must stay the same
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

int **Memory::grow_2d_int_array(int **array, int n1, int n2, const char *name)

{
  if (n1 == 0 || n2 == 0) {
    destroy_2d_int_array(array);
    return NULL;
  }

  if (array == NULL) return create_2d_int_array(n1,n2,name);

  int *data = (int *) srealloc(array[0],n1*n2*sizeof(int),name);
  array = (int **) srealloc(array,n1*sizeof(int *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   create a 2d double array with 2nd index from n2lo to n2hi inclusive 
------------------------------------------------------------------------- */

double **Memory::create_2d_double_array(int n1, int n2lo, int n2hi,
					const char *name)
{
  int n2 = n2hi - n2lo + 1;
  double **array = create_2d_double_array(n1,n2,name);

  for (int i = 0; i < n1; i++) array[i] -= n2lo;
  return array;
}

/* ----------------------------------------------------------------------
   free a 2d double array with 2nd index offset 
------------------------------------------------------------------------- */

void Memory::destroy_2d_double_array(double **array, int offset)
{
  if (array == NULL) return;
  sfree(&array[0][offset]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 3d double array 
------------------------------------------------------------------------- */

double ***Memory::create_3d_double_array(int n1, int n2, int n3,
					 const char *name)
{
  int i,j;

  double *data = (double *) smalloc(n1*n2*n3*sizeof(double),name);
  double **plane = (double **) smalloc(n1*n2*sizeof(double *),name);
  double ***array = (double ***) smalloc(n1*sizeof(double **),name);

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
   free a 3d double array 
------------------------------------------------------------------------- */

void Memory::destroy_3d_double_array(double ***array)
{
  if (array == NULL) return;
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 3d double array
   last 2 dims must stay the same
   if any dim is 0, return NULL 
------------------------------------------------------------------------- */

double ***Memory::grow_3d_double_array(double ***array,
				       int n1, int n2, int n3,
				       const char *name)
{
  int i,j;

  if (n1 == 0 || n2 == 0 || n3 == 0) {
    destroy_3d_double_array(array);
    return NULL;
  }

  if (array == NULL) return create_3d_double_array(n1,n2,n3,name);

  double *data = (double *) srealloc(array[0][0],n1*n2*n3*sizeof(double),name);
  double **plane = (double **) srealloc(array[0],n1*n2*sizeof(double *),name);
  array = (double ***) srealloc(array,n1*sizeof(double **),name);

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
   create a 3d double array with 1st index from n1lo to n1hi inclusive 
------------------------------------------------------------------------- */

double ***Memory::create_3d_double_array(int n1lo, int n1hi, 
					 int n2, int n3, const char *name)
{
  int n1 = n1hi - n1lo + 1;
  double ***array = create_3d_double_array(n1,n2,n3,name);
  return array-n1lo;
}

/* ----------------------------------------------------------------------
   free a 3d double array with 1st index offset 
------------------------------------------------------------------------- */

void Memory::destroy_3d_double_array(double ***array, int offset)
{
  if (array) destroy_3d_double_array(array + offset);
}

/* ----------------------------------------------------------------------
   create a 3d double array with
   1st index from n1lo to n1hi inclusive,
   2nd index from n2lo to n2hi inclusive,
   3rd index from n3lo to n3hi inclusive 
------------------------------------------------------------------------- */

double ***Memory::create_3d_double_array(int n1lo, int n1hi,
					 int n2lo, int n2hi,
					 int n3lo, int n3hi, const char *name)
{
  int n1 = n1hi - n1lo + 1;
  int n2 = n2hi - n2lo + 1;
  int n3 = n3hi - n3lo + 1;
  double ***array = create_3d_double_array(n1,n2,n3,name);

  for (int i = 0; i < n1*n2; i++) array[0][i] -= n3lo;
  for (int i = 0; i < n1; i++) array[i] -= n2lo;
  return array-n1lo;
}

/* ----------------------------------------------------------------------
   free a 3d double array with all 3 indices offset 
------------------------------------------------------------------------- */

void Memory::destroy_3d_double_array(double ***array, int n1_offset,
				     int n2_offset, int n3_offset)
{
  if (array == NULL) return;
  sfree(&array[n1_offset][n2_offset][n3_offset]);
  sfree(&array[n1_offset][n2_offset]);
  sfree(array + n1_offset);
}


/* ----------------------------------------------------------------------
   create a 3d int array 
------------------------------------------------------------------------- */

int ***Memory::create_3d_int_array(int n1, int n2, int n3, const char *name)
{
  int i,j;

  int *data = (int *) smalloc(n1*n2*n3*sizeof(int),name);
  int **plane = (int **) smalloc(n1*n2*sizeof(int *),name);
  int ***array = (int ***) smalloc(n1*sizeof(int **),name);

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
   free a 3d int array 
------------------------------------------------------------------------- */

void Memory::destroy_3d_int_array(int ***array)
{
  if (array == NULL) return;
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 4d double array 
------------------------------------------------------------------------- */

double ****Memory::create_4d_double_array(int n1, int n2, int n3, int n4,
					  const char *name)
{
  int i,j,k;

  double *data = (double *) smalloc(n1*n2*n3*n4*sizeof(double),name);
  double **cube = (double **) smalloc(n1*n2*n3*sizeof(double *),name);
  double ***plane = (double ***) smalloc(n1*n2*sizeof(double **),name);
  double ****array = (double ****) smalloc(n1*sizeof(double ***),name);

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
   free a 4d double array 
------------------------------------------------------------------------- */

void Memory::destroy_4d_double_array(double ****array)
{
  if (array == NULL) return;
  sfree(array[0][0][0]);
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}
