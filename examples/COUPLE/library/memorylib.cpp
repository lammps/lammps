#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include "memorylib.h"
#include "errorlib.h"

/* ---------------------------------------------------------------------- */

MemoryLib::MemoryLib(MPI_Comm comm)
{
  error = new ErrorLib(comm);
}

/* ---------------------------------------------------------------------- */

MemoryLib::~MemoryLib()
{
  delete error;
}

/* ----------------------------------------------------------------------
   safe malloc 
------------------------------------------------------------------------- */

void *MemoryLib::smalloc(int n, const char *name)
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

void MemoryLib::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   safe realloc 
------------------------------------------------------------------------- */

void *MemoryLib::srealloc(void *ptr, int n, const char *name)
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
   create a 2d double array 
------------------------------------------------------------------------- */

double **MemoryLib::create_2d_double_array(int n1, int n2, const char *name)

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
   grow or shrink 1st dim of a 2d double array
   last dim must stay the same
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */

double **MemoryLib::grow_2d_double_array(double **array,
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
   free a 2d double array 
------------------------------------------------------------------------- */

void MemoryLib::destroy_2d_double_array(double **array)

{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}
