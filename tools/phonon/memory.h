#ifndef LMP_MEMORY_H
#define LMP_MEMORY_H

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#include "stdio.h"
#include "stdlib.h"
#include "limits.h"
#include "stdint.h"
#include "inttypes.h"

typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#define ATOBIGINT atoll

class Memory {
 public:
  Memory(){};

  void *smalloc(bigint n, const char *);
  void *srealloc(void *, bigint n, const char *);
  void sfree(void *);
  void fail(const char *);

/* ----------------------------------------------------------------------
   create a 1d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name)
    {
      bigint nbytes = sizeof(TYPE) * n;
      array = (TYPE *) smalloc(nbytes,name);
      return array;
    };

  template <typename TYPE>
    TYPE **create(TYPE **&array, int n, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   grow or shrink 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *grow(TYPE *&array, int n, const char *name) 
    {
      if (array == NULL) return create(array,n,name);
      
      bigint nbytes = sizeof(TYPE) * n;
      array = (TYPE *) srealloc(array,nbytes,name);
      return array;
    };

  template <typename TYPE>
    TYPE **grow(TYPE **&array, int n, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 1d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE *array) 
    {
      sfree(array);
    };

/* ----------------------------------------------------------------------
   create a 1d array with index from nlo to nhi inclusive 
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *create1d_offset(TYPE *&array, int nlo, int nhi, const char *name) 
    {
      bigint nbytes = sizeof(TYPE) * (nhi-nlo+1);
      array = (TYPE *) smalloc(nbytes,name);
      array -= nlo;
      return array;
    }

  template <typename TYPE>
    TYPE **create1d_offset(TYPE **&array, int nlo, int nhi, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 1d array with index offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy1d_offset(TYPE *array, int offset) 
    {
      if (array) sfree(&array[offset]);
    }

/* ----------------------------------------------------------------------
   create a 2d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name) 
    {
      bigint nbytes = sizeof(TYPE) * n1*n2;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = sizeof(TYPE *) * n1;
      array = (TYPE **) smalloc(nbytes,name);
      
      int n = 0;
      for (int i = 0; i < n1; i++) {
	array[i] = &data[n];
	n += n2;
      }
      return array;
    }

  template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **grow(TYPE **&array, int n1, int n2, const char *name) 
    {
      if (array == NULL) return create(array,n1,n2,name);
      
      bigint nbytes = sizeof(TYPE) * n1*n2;
      TYPE *data = (TYPE *) srealloc(array[0],nbytes,name);
      nbytes = sizeof(TYPE *) * n1;
      array = (TYPE **) srealloc(array,nbytes,name);
      
      int n = 0;
      for (int i = 0; i < n1; i++) {
	array[i] = &data[n];
	n += n2;
      }
      return array;
    }

  template <typename TYPE>
    TYPE ***grow(TYPE ***&array, int n1, int n2, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 2d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE **array)
    {
      if (array == NULL) return;
      sfree(array[0]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 2d array with 2nd index from n2lo to n2hi inclusive 
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **create2d_offset(TYPE **&array, int n1, int n2lo, int n2hi,
			   const char *name)
    {
      int n2 = n2hi - n2lo + 1;
      create(array,n1,n2,name);
      for (int i = 0; i < n1; i++) array[i] -= n2lo;
      return array;
    }

  template <typename TYPE>
    TYPE ***create2d_offset(TYPE ***&array, int n1, int n2lo, int n2hi,
			    const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 2d array with 2nd index offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy2d_offset(TYPE **array, int offset)
    {
      if (array == NULL) return;
      sfree(&array[0][offset]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 3d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, int n3, const char *name) 
    {
      bigint nbytes = sizeof(TYPE) * n1*n2*n3;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = sizeof(TYPE *) * n1*n2;
      TYPE **plane = (TYPE **) smalloc(nbytes,name);
      nbytes = sizeof(TYPE **) * n1;
      array = (TYPE ***) smalloc(nbytes,name);
      
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

  template <typename TYPE>
    TYPE ****create(TYPE ****&array, int n1, int n2, int n3, const char *name) 
    {fail(name);}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 3d array
   last 2 dims must stay the same
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***grow(TYPE ***&array, int n1, int n2, int n3, const char *name) 
    {
      if (array == NULL) return create(array,n1,n2,n3,name);
      
      bigint nbytes = sizeof(TYPE) * n1*n2*n3;
      TYPE *data = (TYPE *) srealloc(array[0][0],nbytes,name);
      nbytes = sizeof(TYPE *) * n1*n2;
      TYPE **plane = (TYPE **) srealloc(array[0],nbytes,name);
      nbytes = sizeof(TYPE **) * n1;
      array = (TYPE ***) srealloc(array,nbytes,name);
      
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

  template <typename TYPE>
    TYPE ****grow(TYPE ****&array, int n1, int n2, int n3, const char *name) 
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 3d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ***array) 
    {
      if (array == NULL) return;
      sfree(array[0][0]);
      sfree(array[0]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 3d array with 1st index from n1lo to n1hi inclusive 
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi, 
			    int n2, int n3, const char *name)
    {
      int n1 = n1hi - n1lo + 1;
      create(array,n1,n2,n3,name);
      array -= n1lo;
      return array;
    }

  template <typename TYPE>
    TYPE ****create3d_offset(TYPE ****&array, int n1lo, int n1hi, 
			     int n2, int n3, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   free a 3d array with 1st index offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy3d_offset(TYPE ***array, int offset)
    {
      if (array) destroy(&array[offset]);
    }

/* ----------------------------------------------------------------------
   create a 3d array with
   1st index from n1lo to n1hi inclusive,
   2nd index from n2lo to n2hi inclusive,
   3rd index from n3lo to n3hi inclusive 
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi, 
			    int n2lo, int n2hi, int n3lo, int n3hi,
			    const char *name)
    {
      int n1 = n1hi - n1lo + 1;
      int n2 = n2hi - n2lo + 1;
      int n3 = n3hi - n3lo + 1;
      create(array,n1,n2,n3,name);
      
      for (int i = 0; i < n1*n2; i++) array[0][i] -= n3lo;
      for (int i = 0; i < n1; i++) array[i] -= n2lo;
      array -= n1lo;
      return array;
    }

  template <typename TYPE>
    TYPE ****create3d_offset(TYPE ****&array, int n1lo, int n1hi, 
			     int n2lo, int n2hi, int n3lo, int n3hi,
			     const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   free a 3d array with all 3 indices offset 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy3d_offset(TYPE ***array, 
			  int n1_offset, int n2_offset, int n3_offset)
    {
      if (array == NULL) return;
      sfree(&array[n1_offset][n2_offset][n3_offset]);
      sfree(&array[n1_offset][n2_offset]);
      sfree(&array[n1_offset]);
    }

/* ----------------------------------------------------------------------
   create a 4d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ****create(TYPE ****&array, int n1, int n2, int n3, int n4,
		    const char *name)
    {
      bigint nbytes = sizeof(TYPE) * n1*n2*n3*n4;
      TYPE *data = (double *) smalloc(nbytes,name);
      nbytes = sizeof(TYPE *) * n1*n2*n3;
      TYPE **cube = (double **) smalloc(nbytes,name);
      nbytes = sizeof(TYPE **) * n1*n2;
      TYPE ***plane = (double ***) smalloc(nbytes,name);
      nbytes = sizeof(TYPE ***) * n1;
      array = (double ****) smalloc(nbytes,name);
      
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

  template <typename TYPE>
    TYPE *****create(TYPE *****&array, int n1, int n2, int n3, int n4,
		     const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 4d array 
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ****array)
    {
      if (array == NULL) return;
      sfree(array[0][0][0]);
      sfree(array[0][0]);
      sfree(array[0]);
      sfree(array);
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

#endif
