// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision$
// $Date$
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt
// in the root directory of this source distribution.
// ------------------------------------------------------------- 
#ifndef   __RADIXSORT_H__
#define   __RADIXSORT_H__

#define SORT_CTA_SIZE 256 //This CTA_SIZE must equal 16 * number of radices 

#include "cudpp_globals.h"
#include "cudpp.h"
#include "cudpp_plan.h"


extern "C"
void allocRadixSortStorage(CUDPPRadixSortPlan* plan);

extern "C"
void freeRadixSortStorage(CUDPPRadixSortPlan* plan);

extern "C"
void cudppRadixSortDispatch(void    *keys,
                            void    *values,
                            size_t  numElements,
                            int     keyBits,
                            const   CUDPPRadixSortPlan *plan);


#endif // __RADIXSORT_H__
