// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
//  $Revision: 5633 $
//  $Date: 2009-07-01 15:02:51 +1000 (Wed, 01 Jul 2009) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt 
// in the root directory of this source distribution.
// ------------------------------------------------------------- 

/**
 * @file
 * scan_kernel.cu
 *
 * @brief CUDPP kernel-level scan routines
 */

/** \defgroup cudpp_kernel CUDPP Kernel-Level API
  * The CUDPP Kernel-Level API contains functions that run on the GPU 
  * device across a grid of Cooperative Thread Array (CTA, aka Thread
  * Block).  These kernels are declared \c __global__ so that they 
  * must be invoked from host (CPU) code.  They generally invoke GPU 
  * \c __device__ routines in the CUDPP \link cudpp_cta CTA-Level API\endlink. 
  * Kernel-Level API functions are used by CUDPP 
  * \link cudpp_app Application-Level\endlink functions to implement their 
  * functionality.
  * @{
  */

/** @name Scan Functions
* @{
*/

#include <cudpp_globals.h>
#include "cta/scan_cta.cu"
#include "sharedmem.h"

/**
  * @brief Main scan kernel
  *
  * This __global__ device function performs one level of a multiblock scan on 
  * an arbitrary-dimensioned array in \a d_in, returning the result in \a d_out 
  * (which may point to the same array).  The same function may be used for
  * single or multi-row scans.  To perform a multirow scan, pass the width of 
  * each row of the input row (in elements) in \a dataRowPitch, and the width of 
  * the rows of \a d_blockSums (in elements) in \a blockSumRowPitch, and invoke
  * with a thread block grid with height greater than 1.  
  * 
  * This function peforms one level of a recursive, multiblock scan.  At the 
  * app level, this function is called by cudppScan and cudppMultiScan and used 
  * in combination with vectorAddUniform4() to produce a complete scan.
  *
  * Template parameter \a T is the datatype of the array to be scanned. 
  * Template parameter \a traits is the ScanTraits struct containing 
  * compile-time options for the scan, such as whether it is forward or 
  * backward, exclusive or inclusive, multi- or single-row, etc.
  * 
  * @param[out] d_out The output (scanned) array
  * @param[in]  d_in The input array to be scanned
  * @param[out] d_blockSums The array of per-block sums
  * @param[in]  numElements The number of elements to scan
  * @param[in]  dataRowPitch The width of each row of \a d_in in elements 
  * (for multi-row scans)
  * @param[in]  blockSumRowPitch The with of each row of \a d_blockSums in elements
  * (for multi-row scans)
  */
template<class T, class traits> 
__global__ void scan4(T            *d_out, 
                      const T      *d_in, 
                      T            *d_blockSums, 
                      int          numElements, 
                      unsigned int dataRowPitch,
                      unsigned int blockSumRowPitch)
{
    SharedMemory<T> smem;
    T* temp = smem.getPointer();

    int devOffset, ai, bi, aiDev, biDev;
    T threadScan0[4], threadScan1[4];

    unsigned int blockN = numElements;
    unsigned int blockSumIndex = blockIdx.x;

    if (traits::isMultiRow())
    {
        //int width = __mul24(gridDim.x, blockDim.x) << 1;
        int yIndex     = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;
        devOffset      = __umul24(dataRowPitch, yIndex);
        blockN        += (devOffset << 2);
        devOffset     += __umul24(blockIdx.x, blockDim.x << 1);
        blockSumIndex += __umul24(blockSumRowPitch << 2, yIndex) ;
    }
    else
    {
        devOffset = __umul24(blockIdx.x, (blockDim.x << 1));
    }
    
    // load data into shared memory
    loadSharedChunkFromMem4<T, traits>
        (temp, threadScan0, threadScan1, d_in,
         blockN, devOffset, ai, bi, aiDev, biDev);

    scanCTA<T, traits>(temp, d_blockSums, blockSumIndex);
    
    // write results to device memory
    storeSharedChunkToMem4<T, traits>
        (d_out, threadScan0, threadScan1, temp, 
         blockN, devOffset, ai, bi, aiDev, biDev);

}

/** @} */ // end scan functions
/** @} */ // end cudpp_kernel
