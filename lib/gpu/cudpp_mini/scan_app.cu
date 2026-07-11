// -------------------------------------------------------------
// CUDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 5633 $
// $Date: 2009-07-01 15:02:51 +1000 (Wed, 01 Jul 2009) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt 
// in the root directory of this source distribution.
// ------------------------------------------------------------- 

/**
 * @file
 * scan_app.cu
 *
 * @brief CUDPP application-level scan routines
 */

/** \defgroup cudpp_app CUDPP Application-Level API
  * The CUDPP Application-Level API contains functions
  * that run on the host CPU and invoke GPU routines in 
  * the CUDPP \link cudpp_kernel Kernel-Level API\endlink. 
  * Application-Level API functions are used by
  * CUDPP \link publicInterface Public Interface\endlink
  * functions to implement CUDPP's core functionality.
  * @{
  */

/** @name Scan Functions
 * @{
 */

#include "cudpp.h"
#include "cudpp_util.h"
#include "cudpp_plan.h"
#include "kernel/scan_kernel.cu"
#include "kernel/vector_kernel.cu"


#include <cutil.h>
#include <cstdlib>
#include <cstdio>
#include <assert.h>

/** @brief Perform recursive scan on arbitrary size arrays
  *
  * This is the CPU-side workhorse function of the scan engine.  This function
  * invokes the CUDA kernels which perform the scan on individual blocks. 
  *
  * Scans of large arrays must be split (possibly recursively) into a hierarchy of block scans,
  * where each block is scanned by a single CUDA thread block.  At each recursive level of the
  * scanArrayRecursive first invokes a kernel to scan all blocks of that level, and if the level
  * has more than one block, it calls itself recursively.  On returning from each recursive level,
  * the total sum of each block from the level below is added to all elements of the corresponding
  * block in this level.  See "Parallel Prefix Sum (Scan) in CUDA" for more information (see
  * \ref references ).
  * 
  * Template parameter \a T is the datatype; \a isBackward specifies backward or forward scan; 
  * \a isExclusive specifies exclusive or inclusive scan, and \a op specifies the binary associative
  * operator to be used.
  *
  * @param[out] d_out       The output array for the scan results
  * @param[in]  d_in        The input array to be scanned
  * @param[out] d_blockSums Array of arrays of per-block sums (one array per recursive level, allocated
  *                         by allocScanStorage())
  * @param[in]  numElements The number of elements in the array to scan
  * @param[in]  numRows The number of rows in the array to scan
  * @param[in]  rowPitches  Array of row pitches (one array per recursive level, allocated by 
  *                         allocScanStorage())
  * @param[in]  level       The current recursive level of the scan
  */
template <class T, bool isBackward, bool isExclusive, CUDPPOperator op>
void scanArrayRecursive(T                   *d_out, 
                        const T             *d_in, 
                        T                   **d_blockSums,
                        size_t              numElements,
                        size_t              numRows,
                        const size_t        *rowPitches,
                        int                 level)
{
    unsigned int numBlocks = 
        max(1, (unsigned int)ceil((double)numElements / ((double)SCAN_ELTS_PER_THREAD * CTA_SIZE)));

    unsigned int sharedEltsPerBlock = CTA_SIZE * 2;
      
    unsigned int sharedMemSize = sizeof(T) * sharedEltsPerBlock;

    // divide pitch by four since scan's load/store addresses are for vec4 elements
    unsigned int rowPitch = 1;
    unsigned int blockSumRowPitch = 1;

    if (numRows > 1)
    {
        rowPitch         = rowPitches[level] / 4; 
        blockSumRowPitch = (numBlocks > 1) ? rowPitches[level+1] / 4 : 0;
    }

    bool fullBlock = (numElements == numBlocks * SCAN_ELTS_PER_THREAD * CTA_SIZE);

    // setup execution parameters
    dim3  grid(numBlocks, numRows, 1); 
    dim3  threads(CTA_SIZE, 1, 1);

    // make sure there are no CUDA errors before we start
    CUT_CHECK_ERROR("scanArray before kernels");

    unsigned int traitsCode = 0;
    if (numBlocks > 1) traitsCode |= 1;
    if (numRows > 1)   traitsCode |= 2;
    if (fullBlock)     traitsCode |= 4;

    switch (traitsCode)
    {
    case 0: // single block, single row, non-full block
        scan4<T, ScanTraits<T, op, isBackward, isExclusive, false, false, false> >
               <<< grid, threads, sharedMemSize >>>
               (d_out, d_in, 0, numElements, rowPitch, blockSumRowPitch);
        break;
    case 1: // multiblock, single row, non-full block
        scan4< T, ScanTraits<T, op, isBackward, isExclusive, false, true, false> >
               <<< grid, threads, sharedMemSize >>>
               (d_out, d_in, d_blockSums[level], numElements, rowPitch, blockSumRowPitch);
        break;
    case 2: // single block, multirow, non-full block
        scan4<T, ScanTraits<T, op, isBackward, isExclusive, true, false, false> >
                <<< grid, threads, sharedMemSize >>>
                (d_out, d_in, 0, numElements, rowPitch, blockSumRowPitch);
        break;
    case 3: // multiblock, multirow, non-full block
        scan4<T, ScanTraits<T, op, isBackward, isExclusive, true, true, false> >
                <<< grid, threads, sharedMemSize >>>
                (d_out, d_in, d_blockSums[level], numElements, rowPitch, blockSumRowPitch);
        break;
    case 4: // single block, single row, full block
        scan4<T, ScanTraits<T, op, isBackward, isExclusive, false, false, true> >
               <<< grid, threads, sharedMemSize >>>
               (d_out, d_in, 0, numElements, rowPitch, blockSumRowPitch);
        break;
    case 5: // multiblock, single row, full block
        scan4< T, ScanTraits<T, op, isBackward, isExclusive, false, true, true> >
               <<< grid, threads, sharedMemSize >>>
               (d_out, d_in, d_blockSums[level], numElements, rowPitch, blockSumRowPitch);
        break;
    case 6: // single block, multirow, full block
        scan4<T, ScanTraits<T, op, isBackward, isExclusive, true, false, true> >
                <<< grid, threads, sharedMemSize >>>
                (d_out, d_in, 0, numElements, rowPitch, blockSumRowPitch);
        break;
    case 7: // multiblock, multirow, full block
        scan4<T, ScanTraits<T, op, isBackward, isExclusive, true, true, true> >
                <<< grid, threads, sharedMemSize >>>
                (d_out, d_in, d_blockSums[level], numElements, rowPitch, blockSumRowPitch);
        break;
    }

    CUT_CHECK_ERROR("prescan");

    if (numBlocks > 1)
    {
        // After scanning all the sub-blocks, we are mostly done. But
        // now we need to take all of the last values of the
        // sub-blocks and scan those. This will give us a new value
        // that must be sdded to each block to get the final results.

        scanArrayRecursive<T, isBackward, true, op>
            ((T*)d_blockSums[level], (const T*)d_blockSums[level],
             (T**)d_blockSums, numBlocks, numRows, rowPitches, level + 1); // recursive (CPU) call
        
        vectorAddUniform4<T, op, SCAN_ELTS_PER_THREAD>
            <<< grid, threads >>>(d_out, 
                                  (T*)d_blockSums[level], 
                                  numElements,
                                  rowPitch*4,
                                  blockSumRowPitch*4,
                                  0, 0);
        CUT_CHECK_ERROR("vectorAddUniform");
    }
}

// global
    
#ifdef __cplusplus
extern "C" 
{
#endif

/** @brief Allocate intermediate arrays used by scan.
  *
  * Scans of large arrays must be split (possibly recursively) into a hierarchy 
  * of block scans, where each block is scanned by a single CUDA thread block.  
  * At each recursive level of the scan, we need an array in which to store the 
  * total sums of all blocks in that level.  This function computes the amount 
  * of storage needed and allocates it.
  *
  * @param plan Pointer to CUDPPScanPlan object containing options and number 
  *             of elements, which is used to compute storage requirements, and
  *             within which intermediate storage is allocated.
  */
void allocScanStorage(CUDPPScanPlan *plan)
{
    //assert(config->_numEltsAllocated == 0); // shouldn't be called 

    plan->m_numEltsAllocated = plan->m_numElements;

    size_t numElts = plan->m_numElements;
    
    size_t level = 0;

    do
    {       
        size_t numBlocks = 
            max(1, (unsigned int)ceil((double)numElts / ((double)SCAN_ELTS_PER_THREAD * CTA_SIZE)));
        if (numBlocks > 1)
        {
            level++;
        }
        numElts = numBlocks;
    } while (numElts > 1);

    size_t elementSize = 0;

    switch(plan->m_config.datatype)
    {
    case CUDPP_INT:
        plan->m_blockSums = (void**) malloc(level * sizeof(int*));
        elementSize = sizeof(int);
        break;
    case CUDPP_UINT:
        plan->m_blockSums = (void**) malloc(level * sizeof(unsigned int*));
        elementSize = sizeof(unsigned int);
        break;
    case CUDPP_FLOAT:
        plan->m_blockSums = (void**) malloc(level * sizeof(float*));
        elementSize = sizeof(float);
        break;
    default:
        break;
    }

    plan->m_numLevelsAllocated = level;
    numElts = plan->m_numElements;
    size_t numRows = plan->m_numRows;
    plan->m_numRowsAllocated = numRows;
    plan->m_rowPitches = 0;

    if (numRows > 1)
    {
        plan->m_rowPitches = (size_t*) malloc((level + 1) * sizeof(size_t));
        plan->m_rowPitches[0] = plan->m_rowPitch;
    }

    level = 0;

    do
    {       
        size_t numBlocks = 
            max(1, (unsigned int)ceil((double)numElts / ((double)SCAN_ELTS_PER_THREAD * CTA_SIZE)));
        if (numBlocks > 1) 
        {
            // Use cudaMallocPitch for multi-row block sums to ensure alignment
            if (numRows > 1)
            {
                size_t dpitch;
                CUDA_SAFE_CALL( cudaMallocPitch((void**) &(plan->m_blockSums[level]), 
                                                &dpitch,
                                                numBlocks * elementSize, 
                                                numRows));
                plan->m_rowPitches[level+1] = dpitch / elementSize;
                level++;
            }
            else
            {
                CUDA_SAFE_CALL(cudaMalloc((void**) &(plan->m_blockSums[level++]),  
                                          numBlocks * elementSize));
            }
        }
        numElts = numBlocks;
    } while (numElts > 1);

    CUT_CHECK_ERROR("allocScanStorage");
}

/** @brief Deallocate intermediate block sums arrays in a CUDPPScanPlan object.
  *
  * These arrays must have been allocated by allocScanStorage(), which is called
  * by the constructor of cudppScanPlan().  
  *
  * @param plan Pointer to CUDPPScanPlan object initialized by allocScanStorage().
  */
void freeScanStorage(CUDPPScanPlan *plan)
{
    for (unsigned int i = 0; i < plan->m_numLevelsAllocated; i++)
    {
        cudaFree(plan->m_blockSums[i]);
    }

    CUT_CHECK_ERROR("freeScanStorage");

    free((void**)plan->m_blockSums);
    if (plan->m_numRows > 1)
        free((void*)plan->m_rowPitches);

    plan->m_blockSums = 0;
    plan->m_numEltsAllocated = 0;
    plan->m_numLevelsAllocated = 0;
}


/** @brief Dispatch function to perform a scan (prefix sum) on an
  * array with the specified configuration.
  *
  * This is the dispatch routine which calls scanArrayRecursive() with 
  * appropriate template parameters and arguments to achieve the scan as 
  * specified in \a plan. 
  * 
  * @param[out] d_out    The output array of scan results
  * @param[in]  d_in     The input array
  * @param[in]  numElements The number of elements to scan
  * @param[in]  numRows     The number of rows to scan in parallel
  * @param[in]  plan     Pointer to CUDPPScanPlan object containing scan options
  *                      and intermediate storage
  */
void cudppScanDispatch(void                *d_out, 
                       const void          *d_in, 
                       size_t              numElements,
                       size_t              numRows,
                       const CUDPPScanPlan *plan)
{    
    if (CUDPP_OPTION_EXCLUSIVE & plan->m_config.options)
    {
        if (CUDPP_OPTION_BACKWARD & plan->m_config.options)
        {
            switch (plan->m_config.datatype)
            {
            case CUDPP_INT:

                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<int, true, true, CUDPP_ADD>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<int, true, true, CUDPP_MULTIPLY>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<int, true, true, CUDPP_MAX>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<int, true, true, CUDPP_MIN>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }
              
                break;

            case CUDPP_UINT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:                 
                    scanArrayRecursive<unsigned int, true, true, CUDPP_ADD>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:                 
                    scanArrayRecursive<unsigned int, true, true, CUDPP_MULTIPLY>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<unsigned int, true, true, CUDPP_MAX>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<unsigned int, true, true, CUDPP_MIN>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }

                break;

            case CUDPP_FLOAT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<float, true, true,  CUDPP_ADD>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<float, true, true,  CUDPP_MULTIPLY>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<float, true, true, CUDPP_MAX>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<float, true, true, CUDPP_MIN>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }
                break; 

            default:
                break; 
            }
        }
        else
        {
            switch (plan->m_config.datatype)
            {
            case CUDPP_INT:

                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<int, false, true, CUDPP_ADD>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<int, false, true, CUDPP_MULTIPLY>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<int, false, true, CUDPP_MAX>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<int, false, true, CUDPP_MIN>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }

                break;
                    
            case CUDPP_UINT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:                 
                    scanArrayRecursive<unsigned int, false, true, CUDPP_ADD>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:                 
                    scanArrayRecursive<unsigned int, false, true, CUDPP_MULTIPLY>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<unsigned int, false, true, CUDPP_MAX>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<unsigned int, false, true, CUDPP_MIN>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                            
                }
        
                break;       
            
            case CUDPP_FLOAT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<float, false, true, CUDPP_ADD>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<float, false, true, CUDPP_MULTIPLY>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<float, false, true, CUDPP_MAX>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<float, false, true, CUDPP_MIN>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }            
                break;

            default:
                break; 
            }
        }
    }
    else
    {
        if (CUDPP_OPTION_BACKWARD & plan->m_config.options)
        {
            switch (plan->m_config.datatype)
            {
            case CUDPP_INT:

                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<int, true, false, CUDPP_ADD>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<int, true, false, CUDPP_MULTIPLY>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<int, true, false, CUDPP_MAX>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<int, true, false, CUDPP_MIN>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }
              
                break;

            case CUDPP_UINT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:                 
                    scanArrayRecursive<unsigned int, true, false, CUDPP_ADD>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:                 
                    scanArrayRecursive<unsigned int, true, false, CUDPP_MULTIPLY>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<unsigned int, true, false, CUDPP_MAX>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<unsigned int, true, false, CUDPP_MIN>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }

                break;

            case CUDPP_FLOAT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<float, true, false, CUDPP_ADD>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<float, true, false, CUDPP_MULTIPLY>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<float, true, false, CUDPP_MAX>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<float, true, false, CUDPP_MIN>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }
                break; 

            default:
                break; 
            }
        }
        else
        {
            switch (plan->m_config.datatype)
            {
            case CUDPP_INT:

                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<int, false, false, CUDPP_ADD>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<int, false, false, CUDPP_MULTIPLY>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<int, false, false, CUDPP_MAX>
                        ((int*)d_out, (const int*)d_in, 
                         (int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<int, false, false, CUDPP_MIN>
                        ((int*)d_out, (const int*)d_in, 
                        (int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }

                break;
                    
            case CUDPP_UINT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:                 
                    scanArrayRecursive<unsigned int, false, false, CUDPP_ADD>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:                 
                    scanArrayRecursive<unsigned int, false, false, CUDPP_MULTIPLY>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<unsigned int, false, false, CUDPP_MAX>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                         (unsigned int**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<unsigned int, false, false, CUDPP_MIN>
                        ((unsigned int*)d_out, (const unsigned int*)d_in, 
                        (unsigned int**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                            
                }
        
                break;       
            
            case CUDPP_FLOAT:
                switch(plan->m_config.op)
                {
                case CUDPP_ADD:
                    scanArrayRecursive<float, false, false, CUDPP_ADD>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MULTIPLY:
                    scanArrayRecursive<float, false, false, CUDPP_MULTIPLY>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MAX:
                    scanArrayRecursive<float, false, false, CUDPP_MAX>
                        ((float*)d_out, (const float*)d_in, 
                         (float**)plan->m_blockSums, 
                         numElements, numRows, plan->m_rowPitches, 0);
                    break;
                case CUDPP_MIN:
                    scanArrayRecursive<float, false, false, CUDPP_MIN>
                        ((float*)d_out, (const float*)d_in, 
                        (float**)plan->m_blockSums, 
                        numElements, numRows, plan->m_rowPitches, 0);
                    break;
                default:
                    break;
                }            
                break;

            default:
                break; 
            }
        }  
    }
}

#ifdef __cplusplus
}
#endif

/** @} */ // end scan functions
/** @} */ // end cudpp_app
