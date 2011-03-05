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

void *Memory::smalloc(bigint n, const char *name)
{
  if (n == 0) return NULL;
#if defined(MALLOC_MEMALIGN)
  void *ptr;
  posix_memalign(&ptr, MALLOC_MEMALIGN, n);
#else
  void *ptr = malloc(n);
#endif
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate " BIGINT_FORMAT 
	    " bytes for array %s",n,name);
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

void *Memory::srealloc(void *ptr, bigint n, const char *name)
{
  if (n == 0) {
    sfree(ptr);
    return NULL;
  }

  ptr = realloc(ptr,n);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to reallocate " BIGINT_FORMAT 
	    " bytes for array %s",n,name);
    error->one(str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   create a 1d double array with index from nlo to nhi inclusive 
------------------------------------------------------------------------- */

double *Memory::create_1d_double_array(int nlo, int nhi, const char *name)
{
  double *array;
  if (create_1d_array(&array,nlo,nhi,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 1d double array with index offset 
------------------------------------------------------------------------- */

void Memory::destroy_1d_double_array(double *array, int offset)
{
  destroy_1d_array(array, offset);
}

/* ----------------------------------------------------------------------
   create a 2d double array 
------------------------------------------------------------------------- */

double **Memory::create_2d_double_array(int n1, int n2, const char *name)
{
  double **array;
  if (create_2d_array(&array,n1, n2, name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 2d double array 
------------------------------------------------------------------------- */

void Memory::destroy_2d_double_array(double **array)
{
  destroy_2d_array(array);
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d double array
   last dim must stay the same
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

double **Memory::grow_2d_double_array(double **array,
				      int n1, int n2, const char *name)
{
  if (grow_2d_array(&array, n1, n2, name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   create a 2d int array
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

int **Memory::create_2d_int_array(int n1, int n2, const char *name)
{
  int **array;
  if (create_2d_array(&array,n1,n2,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 2d int array 
------------------------------------------------------------------------- */

void Memory::destroy_2d_int_array(int **array)
{
  destroy_2d_array(array);
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d int array
   last dim must stay the same
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

int **Memory::grow_2d_int_array(int **array, int n1, int n2, const char *name)
{
  if (grow_2d_array(&array, n1, n2, name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   create a 2d double array with 2nd index from n2lo to n2hi inclusive 
------------------------------------------------------------------------- */

double **Memory::create_2d_double_array(int n1, int n2lo, int n2hi,
					const char *name)
{
  double **array;
  if (create_2d_array(&array, n1, n2lo, n2hi, name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 2d double array with 2nd index offset 
------------------------------------------------------------------------- */

void Memory::destroy_2d_double_array(double **array, int offset)
{
  destroy_2d_array(array, offset);
}

/* ----------------------------------------------------------------------
   create a 2d float array 
------------------------------------------------------------------------- */

float **Memory::create_2d_float_array(int n1, int n2, const char *name)
{
  float **array;
  if (create_2d_array(&array, n1, n2, name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 2d float array 
------------------------------------------------------------------------- */

void Memory::destroy_2d_float_array(float **array)
{
  destroy_2d_array(array);
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d float array
   last dim must stay the same
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

float **Memory::grow_2d_float_array(float **array,
				      int n1, int n2, const char *name)

{
  if (grow_2d_array(&array, n1, n2, name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   create a 2d float array with 2nd index from n2lo to n2hi inclusive 
------------------------------------------------------------------------- */

float **Memory::create_2d_float_array(int n1, int n2lo, int n2hi,
					const char *name)
{
  float **array;
  if (create_2d_array(&array, n1, n2lo, n2hi, name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 2d float array with 2nd index offset 
------------------------------------------------------------------------- */

void Memory::destroy_2d_float_array(float **array, int offset)
{
  destroy_2d_array(array,offset);
}

/* ----------------------------------------------------------------------
   create a 3d double array 
------------------------------------------------------------------------- */

double ***Memory::create_3d_double_array(int n1, int n2, int n3,
					 const char *name)
{
  double ***array;
  if (create_3d_array(&array,n1,n2,n3,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 3d double array 
------------------------------------------------------------------------- */

void Memory::destroy_3d_double_array(double ***array)
{
  destroy_3d_array(array);
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
  if (grow_3d_array(&array,n1,n2,n3,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   create a 3d double array with 1st index from n1lo to n1hi inclusive 
------------------------------------------------------------------------- */

double ***Memory::create_3d_double_array(int n1lo, int n1hi, 
					 int n2, int n3, const char *name)
{
  double ***array;
  if (create_3d_array(&array,n1lo,n1hi,n2,n3,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 3d double array with 1st index offset 
------------------------------------------------------------------------- */

void Memory::destroy_3d_double_array(double ***array, int offset)
{
  destroy_3d_array(array,offset);
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
  double ***array;
  if (create_3d_array(&array,n1lo,n1hi,n2lo,n2hi,n3lo,n3hi,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 3d double array with all 3 indices offset 
------------------------------------------------------------------------- */

void Memory::destroy_3d_double_array(double ***array, int n1_offset,
				     int n2_offset, int n3_offset)
{
  destroy_3d_array(array,n1_offset,n2_offset,n3_offset);
}

/* ----------------------------------------------------------------------
   create a 3d float array 
------------------------------------------------------------------------- */

float ***Memory::create_3d_float_array(int n1, int n2, int n3,
					 const char *name)
{
  float ***array;
  if (create_3d_array(&array,n1,n2,n3,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 3d float array 
------------------------------------------------------------------------- */

void Memory::destroy_3d_float_array(float ***array)
{
  destroy_3d_array(array);
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 3d float array
   last 2 dims must stay the same
   if any dim is 0, return NULL 
------------------------------------------------------------------------- */

float ***Memory::grow_3d_float_array(float ***array,
				       int n1, int n2, int n3,
				       const char *name)
{ 
  if (grow_3d_array(&array,n1,n2,n3,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   create a 3d float array with 1st index from n1lo to n1hi inclusive 
------------------------------------------------------------------------- */

float ***Memory::create_3d_float_array(int n1lo, int n1hi, 
					 int n2, int n3, const char *name)
{
  float ***array;
  if (create_3d_array(&array,n1lo,n1hi,n2,n3,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 3d float array with 1st index offset 
------------------------------------------------------------------------- */

void Memory::destroy_3d_float_array(float ***array, int offset)
{
  destroy_3d_array(array,offset);
}

/* ----------------------------------------------------------------------
   create a 3d float array with
   1st index from n1lo to n1hi inclusive,
   2nd index from n2lo to n2hi inclusive,
   3rd index from n3lo to n3hi inclusive 
------------------------------------------------------------------------- */

float ***Memory::create_3d_float_array(int n1lo, int n1hi,
					 int n2lo, int n2hi,
					 int n3lo, int n3hi, const char *name)
{
  float ***array;
  if (create_3d_array(&array,n1lo,n1hi,n2lo,n2hi,n3lo,n3hi,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 3d float array with all 3 indices offset 
------------------------------------------------------------------------- */

void Memory::destroy_3d_float_array(float ***array, int n1_offset,
				     int n2_offset, int n3_offset)
{
  destroy_3d_array(array,n1_offset, n2_offset, n3_offset);
}

/* ----------------------------------------------------------------------
   create a 3d int array 
------------------------------------------------------------------------- */

int ***Memory::create_3d_int_array(int n1, int n2, int n3, const char *name)
{
  int ***array;
  if (create_3d_array(&array,n1,n2,n3,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 3d int array 
------------------------------------------------------------------------- */

void Memory::destroy_3d_int_array(int ***array)
{
  destroy_3d_array(array);
}

/* ----------------------------------------------------------------------
   create a 4d double array 
------------------------------------------------------------------------- */

double ****Memory::create_4d_double_array(int n1, int n2, int n3, int n4,
					  const char *name)
{
  double ****array;
  if (create_4d_array(&array,n1,n2,n3,n4,name))
    return array;
  else return NULL;
}

/* ----------------------------------------------------------------------
   free a 4d double array 
------------------------------------------------------------------------- */

void Memory::destroy_4d_double_array(double ****array)
{
  destroy_4d_array(array);
}
