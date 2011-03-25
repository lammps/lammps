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

#include "lmptype.h"
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

void *Memory::smalloc(bigint nbytes, const char *name)
{
  if (nbytes == 0) return NULL;

  void *ptr = malloc(nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate " BIGINT_FORMAT "bytes for array %s",
	    nbytes,name);
    error->one(str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc 
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
  if (nbytes == 0) {
    sfree(ptr);
    return NULL;
  }

  ptr = realloc(ptr,nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to reallocate " BIGINT_FORMAT "bytes for array %s",
	    nbytes,name);
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
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(const char *name)
{
  char str[128];
  sprintf(str,"Cannot create/grow a vector/array of pointers for %s",name);
  error->one(str);
}

/* ----------------------------------------------------------------------
   older routines, will be deprecated at some point
------------------------------------------------------------------------- */

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
