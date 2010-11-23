// -------------------------------------------------------------
// CUDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision$
// $Date$
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt 
// in the root directory of this source distribution.
// ------------------------------------------------------------- 

/**
 * @file
 * radixsort_app.cu
 *   
 * @brief CUDPP application-level radix sorting routines
 */

/** @addtogroup cudpp_app 
 * @{
 */

/** @name RadixSort Functions
 * @{
 */
 

#include "cudpp.h"
#include "cudpp_util.h"
#include "cudpp_radixsort.h"
#include "cudpp_scan.h"
#include "kernel/radixsort_kernel.cu"

#include <cutil.h>
#include <cstdlib>
#include <cstdio>
#include <assert.h>

typedef unsigned int uint;

/** @brief Perform one step of the radix sort.  Sorts by nbits key bits per step, 
* starting at startbit.
* 
* Uses cudppScanDispatch() for the prefix sum of radix counters.
* 
* @param[in,out] keys Keys to be sorted.
* @param[in,out] values Associated values to be sorted (through keys).
* @param[in] plan Configuration information for RadixSort.
* @param[in] numElements Number of elements in the sort.
**/
template<uint nbits, uint startbit, bool flip, bool unflip>
void radixSortStep(uint *keys, 
                   uint *values, 
                   const CUDPPRadixSortPlan *plan,
                   uint numElements)
{
    const uint eltsPerBlock = SORT_CTA_SIZE * 4;
    const uint eltsPerBlock2 = SORT_CTA_SIZE * 2;

    bool fullBlocks = ((numElements % eltsPerBlock) == 0);
    uint numBlocks = (fullBlocks) ? 
        (numElements / eltsPerBlock) : 
    (numElements / eltsPerBlock + 1);
    uint numBlocks2 = ((numElements % eltsPerBlock2) == 0) ?
        (numElements / eltsPerBlock2) : 
    (numElements / eltsPerBlock2 + 1);

    bool loop = numBlocks > 65535;
    uint blocks = loop ? 65535 : numBlocks;
    uint blocksFind = loop ? 65535 : numBlocks2;
    uint blocksReorder = loop ? 65535 : numBlocks2;

    uint threshold = fullBlocks ? plan->m_persistentCTAThresholdFullBlocks[0] : plan->m_persistentCTAThreshold[0];

    bool persist = plan->m_bUsePersistentCTAs && (numElements >= threshold);

    if (persist)
    {
        loop = (numElements > 262144) || (numElements >= 32768 && numElements < 65536);
        
        blocks = numBlocks;
        blocksFind = numBlocks2;
        blocksReorder = numBlocks2;

        // Run an empty kernel -- this seems to reset some of the CTA scheduling hardware
        // on GT200, resulting in better scheduling and lower run times
        if (startbit > 0)
        {
            emptyKernel<<<numCTAs(emptyKernel), SORT_CTA_SIZE>>>();
        }
    }

    if (fullBlocks)
    {
        if (loop)
        {
            if (persist) 
            {
                blocks = flip? numCTAs(radixSortBlocks<4, 0, true, true, true>) : 
                               numCTAs(radixSortBlocks<4, 0, true, false, true>);
            }

            radixSortBlocks<nbits, startbit, true, flip, true>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)plan->m_tempValues, (uint4*)keys, (uint4*)values, numElements, numBlocks);
        }
        else
        {
            radixSortBlocks<nbits, startbit, true, flip, false>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)plan->m_tempValues, (uint4*)keys, (uint4*)values, numElements, numBlocks);
        }
    }
    else
    {
        if (loop)
        {
            if (persist) 
            {
                blocks = flip ? numCTAs(radixSortBlocks<4, 0, false, true, true>) : 
                                numCTAs(radixSortBlocks<4, 0, false, false, true>);
            }

            radixSortBlocks<nbits, startbit, false, flip, true>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)plan->m_tempValues, (uint4*)keys, (uint4*)values, numElements, numBlocks);
        }
        else
        {
            radixSortBlocks<nbits, startbit, false, flip, false>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)plan->m_tempValues, (uint4*)keys, (uint4*)values, numElements, numBlocks);
        }
    }

    CUT_CHECK_ERROR("radixSortBlocks");

    if (fullBlocks)
    {
        if (loop)
        {
            if (persist) 
            {
                blocksFind = numCTAs(findRadixOffsets<0, true, true>);
            }
            findRadixOffsets<startbit, true, true>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);
        }
        else
        {
            findRadixOffsets<startbit, true, false>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);
        }
    }
    else
    {
        if (loop)
        {
            if (persist) 
            {
                blocksFind = numCTAs(findRadixOffsets<0, false, true>);
            }
            findRadixOffsets<startbit, false, true>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);
        }
        else
        {
            findRadixOffsets<startbit, false, false>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);
        }
    }

    CUT_CHECK_ERROR("findRadixOffsets");

    cudppScanDispatch(plan->m_countersSum, plan->m_counters, 16*numBlocks2, 1, plan->m_scanPlan);

    if (fullBlocks)
    {
        if (plan->m_bManualCoalesce)
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ? numCTAs(reorderData<0, true, true, true, true>) :
                                             numCTAs(reorderData<0, true, true, false, true>);
                }
                reorderData<startbit, true, true, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
            else
            {
                reorderData<startbit, true, true, unflip, false>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
        }
        else
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ? numCTAs(reorderData<0, true, false, true, true>) :
                                             numCTAs(reorderData<0, true, false, false, true>);
                }
                reorderData<startbit, true, false, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
            else
            {
                reorderData<startbit, true, false, unflip, false>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
        }
    }
    else
    {
        if (plan->m_bManualCoalesce)
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ? 
                        numCTAs(reorderData<0, false, true, true, true>) :
                        numCTAs(reorderData<0, false, true, false, true>);
                }
                reorderData<startbit, false, true, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
            else
            {
                reorderData<startbit, false, true, unflip, false>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
        }
        else
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ?
                        numCTAs(reorderData<0, false, false, true, true>) :
                        numCTAs(reorderData<0, false, false, false, true>);
                }
                reorderData<startbit, false, false, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
            else
            {
                reorderData<startbit, false, false, unflip, false>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, values, (uint2*)plan->m_tempKeys, (uint2*)plan->m_tempValues, 
                    plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, numElements, numBlocks2);
            }
        }
    }

    CUT_CHECK_ERROR("radixSortStep");
}

/**
 * @brief Single-block optimization for sorts of fewer than 4 * CTA_SIZE elements
 * 
 * @param[in,out] keys  Keys to be sorted.
 * @param[in,out] values Associated values to be sorted (through keys).
 * @param numElements Number of elements in the sort.
**/
template <bool flip>
void radixSortSingleBlock(uint *keys, 
                          uint *values, 
                          uint numElements)
{
    bool fullBlocks = (numElements % (SORT_CTA_SIZE * 4) == 0);
    if (fullBlocks)
    {
        radixSortBlocks<32, 0, true, flip, false>
            <<<1, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)keys, (uint4*)values, 
                 (uint4*)keys, (uint4*)values, 
                 numElements, 0);
    }
    else
    {
        radixSortBlocks<32, 0, false, flip, false>
            <<<1, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)keys, (uint4*)values, 
                 (uint4*)keys, (uint4*)values, 
                 numElements, 0);
    }

    if (flip) unflipFloats<<<1, SORT_CTA_SIZE>>>(keys, numElements);

    CUT_CHECK_ERROR("radixSortSingleBlock");
}

/**
 * @brief Main radix sort function
 * 
 * Main radix sort function.  Sorts in place in the keys and values arrays,
 * but uses the other device arrays as temporary storage.  All pointer 
 * parameters are device pointers.  Uses cudppScan() for the prefix sum of 
 * radix counters.
 * 
 * @param[in,out] keys Keys to be sorted.
 * @param[in,out] values Associated values to be sorted (through keys).
 * @param[in] plan Configuration information for RadixSort.
 * @param[in] numElements Number of elements in the sort.
 * @param[in] flipBits Is set true if key datatype is a float 
 *                 (neg. numbers) for special float sorting operations.
 * @param[in] keyBits Number of interesting bits in the key
 **/
void radixSort(uint *keys,                         
               uint* values,               
               const CUDPPRadixSortPlan *plan,               
               size_t numElements,
               bool flipBits,
               int keyBits)
{
    if(numElements <= WARP_SIZE)
    {
        if (flipBits)
            radixSortSingleWarp<true><<<1, numElements>>>
                (keys, values, numElements);
        else
            radixSortSingleWarp<false><<<1, numElements>>>
                (keys, values, numElements);

        CUT_CHECK_ERROR("radixSortSingleWarp");        
        return;
    }
#ifdef __DEVICE_EMULATION__
    printf("bits: %d\n", keyBits);
#endif
    
    if(numElements <= SORT_CTA_SIZE * 4)
    {
        if (flipBits)
            radixSortSingleBlock<true>(keys, values, numElements);
        else
            radixSortSingleBlock<false>(keys, values, numElements);
        return;
    }
        
    // flip float bits on the first pass, unflip on the last pass    
    if (flipBits) 
    {               
        radixSortStep<4,  0, true, false>
            (keys, values, plan, numElements);            
    }
    else
    {     
        radixSortStep<4,  0, false, false>
            (keys, values, plan, numElements);           
    }

    if (keyBits > 4)
    {                   
        radixSortStep<4,  4, false, false>
            (keys, values, plan, numElements);            
    }
    if (keyBits > 8)
    {                                   
        radixSortStep<4,  8, false, false>
            (keys, values, plan, numElements);            
    }
    if (keyBits > 12)
    {                   
        radixSortStep<4, 12, false, false>
            (keys, values, plan, numElements);            
    }
    if (keyBits > 16)
    {                   
        radixSortStep<4, 16, false, false>
            (keys, values, plan, numElements);            
    }
    if (keyBits > 20)
    {                   
        radixSortStep<4, 20, false, false>
            (keys, values, plan, numElements);            
    }
    if (keyBits > 24)
    {                   
        radixSortStep<4, 24, false, false>
            (keys, values, plan, numElements);         
    }
    if (keyBits > 28)
    {
        if (flipBits) // last pass
        {                       
            radixSortStep<4, 28, false, true>
                (keys, values, plan, numElements);
        }
        else
        {                       
            radixSortStep<4, 28, false, false>
                (keys, values, plan, numElements);            
        }
    }
}

/**
 * @brief Wrapper to call main radix sort function. For float configuration.
 * 
 * Calls the main radix sort function. For float configuration.
 * 
 * @param[in,out] keys Keys to be sorted.
 * @param[in,out] values Associated values to be sorted (through keys).
 * @param[in] plan Configuration information for RadixSort.
 * @param[in] numElements Number of elements in the sort.
 * @param[in] negativeKeys Is set true if key datatype has neg. numbers.
 * @param[in] keyBits Number of interesting bits in the key
 **/
extern "C"
void radixSortFloatKeys(float* keys, 
                        uint* values, 
                        const CUDPPRadixSortPlan *plan,
                        size_t numElements,            
                        bool  negativeKeys,
                        int keyBits)
{
   
    radixSort((uint*)keys, (uint*)values, plan, 
              numElements, negativeKeys, keyBits);
}

/** @brief Perform one step of the radix sort.  Sorts by nbits key bits per step, 
 * starting at startbit.
 * 
 * @param[in,out] keys  Keys to be sorted.
 * @param[in] plan Configuration information for RadixSort.
 * @param[in] numElements Number of elements in the sort. 
**/
template<uint nbits, uint startbit, bool flip, bool unflip>
void radixSortStepKeysOnly(uint *keys, 
                           const CUDPPRadixSortPlan *plan,
                           uint numElements)
{
    const uint eltsPerBlock = SORT_CTA_SIZE * 4;
    const uint eltsPerBlock2 = SORT_CTA_SIZE * 2;

    bool fullBlocks = ((numElements % eltsPerBlock) == 0);
    uint numBlocks = (fullBlocks) ? 
        (numElements / eltsPerBlock) : 
    (numElements / eltsPerBlock + 1);
    uint numBlocks2 = ((numElements % eltsPerBlock2) == 0) ?
        (numElements / eltsPerBlock2) : 
    (numElements / eltsPerBlock2 + 1);

    bool loop = numBlocks > 65535;
    
    uint blocks = loop ? 65535 : numBlocks;
    uint blocksFind = loop ? 65535 : numBlocks2;
    uint blocksReorder = loop ? 65535 : numBlocks2;

    uint threshold = fullBlocks ? plan->m_persistentCTAThresholdFullBlocks[1] : plan->m_persistentCTAThreshold[1];

    bool persist = plan->m_bUsePersistentCTAs && (numElements >= threshold);

    if (persist)
    {
        loop = (numElements > 262144) || (numElements >= 32768 && numElements < 65536);
        
        blocks = numBlocks;
        blocksFind = numBlocks2;
        blocksReorder = numBlocks2;
    }

    if (fullBlocks)
    {
        if (loop)
        {
            if (persist) 
            {
                blocks = flip ? numCTAs(radixSortBlocksKeysOnly<4, 0, true, true, true>) : 
                                numCTAs(radixSortBlocksKeysOnly<4, 0, true, false, true>);
            }

            radixSortBlocksKeysOnly<nbits, startbit, true, flip, true>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)keys, numElements, numBlocks);
        }
        else
            radixSortBlocksKeysOnly<nbits, startbit, true, flip, false>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)keys, numElements, numBlocks);
    }
    else
    {
        if (loop)
        {
            if (persist) 
            {
                blocks = flip ? numCTAs(radixSortBlocksKeysOnly<4, 0, false, true, true>) : 
                                numCTAs(radixSortBlocksKeysOnly<4, 0, false, false, true>);
            }

            radixSortBlocksKeysOnly<nbits, startbit, false, flip, true>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)keys, numElements, numBlocks);
        }
        else
            radixSortBlocksKeysOnly<nbits, startbit, false, flip, false>
                <<<blocks, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint4*)plan->m_tempKeys, (uint4*)keys, numElements, numBlocks);

    }

    if (fullBlocks)
    {
        if (loop)
        {
            if (persist) 
            {
                blocksFind = numCTAs(findRadixOffsets<0, true, true>);
            }
            findRadixOffsets<startbit, true, true>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);
        }
        else
            findRadixOffsets<startbit, true, false>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);
    }
    else
    {
        if (loop)
        {
            if (persist) 
            {
                blocksFind = numCTAs(findRadixOffsets<0, false, true>);
            }
            findRadixOffsets<startbit, false, true>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);
        }
        else
            findRadixOffsets<startbit, false, false>
                <<<blocksFind, SORT_CTA_SIZE, 3 * SORT_CTA_SIZE * sizeof(uint)>>>
                ((uint2*)plan->m_tempKeys, plan->m_counters, plan->m_blockOffsets, numElements, numBlocks2);

    }

    cudppScanDispatch(plan->m_countersSum, plan->m_counters, 16*numBlocks2, 1, plan->m_scanPlan);

    if (fullBlocks)
    {
        if (plan->m_bManualCoalesce)
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ? 
                        numCTAs(reorderDataKeysOnly<0, true, true, true, true>) : 
                        numCTAs(reorderDataKeysOnly<0, true, true, false, true>);
                }
                reorderDataKeysOnly<startbit, true, true, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                    numElements, numBlocks2);
            }
            else
                reorderDataKeysOnly<startbit, true, true, unflip, false>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                     numElements, numBlocks2);
        }
        else
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ?
                        numCTAs(reorderDataKeysOnly<0, true, false, true, true>) :
                        numCTAs(reorderDataKeysOnly<0, true, false, false, true>);
                }
                reorderDataKeysOnly<startbit, true, false, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                    numElements, numBlocks2);
            }
            else
                reorderDataKeysOnly<startbit, true, false, unflip, false>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                     numElements, numBlocks2);
        }
    }
    else
    {
        if (plan->m_bManualCoalesce)
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ? 
                        numCTAs(reorderDataKeysOnly<0, false, true, true, true>) :
                        numCTAs(reorderDataKeysOnly<0, false, true, false, true>);
                }
                reorderDataKeysOnly<startbit, false, true, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                    numElements, numBlocks2);
            }
            else
                reorderDataKeysOnly<startbit, false, true, unflip, false>
                <<<blocksReorder, SORT_CTA_SIZE>>>
                (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                numElements, numBlocks2);
        }
        else
        {
            if (loop)
            {
                if (persist) 
                {
                    blocksReorder = unflip ?
                        numCTAs(reorderDataKeysOnly<0, false, false, true, true>) :
                        numCTAs(reorderDataKeysOnly<0, false, false, false, true>);
                }
                reorderDataKeysOnly<startbit, false, false, unflip, true>
                    <<<blocksReorder, SORT_CTA_SIZE>>>
                    (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                    numElements, numBlocks2);
            }
            else
                reorderDataKeysOnly<startbit, false, false, unflip, false>
                <<<blocksReorder, SORT_CTA_SIZE>>>
                (keys, (uint2*)plan->m_tempKeys, plan->m_blockOffsets, plan->m_countersSum, plan->m_counters, 
                numElements, numBlocks2);
        }
    }

    CUT_CHECK_ERROR("radixSortStepKeysOnly");
}

/**
 * @brief Optimization for sorts of fewer than 4 * CTA_SIZE elements (keys only).
 * 
 * @param[in,out] keys Keys to be sorted.
 * @param numElements Number of elements in the sort.
**/
template <bool flip>
void radixSortSingleBlockKeysOnly(uint *keys, 
                                  uint numElements)
{
    bool fullBlocks = (numElements % (SORT_CTA_SIZE * 4) == 0);
    if (fullBlocks)
    {
        radixSortBlocksKeysOnly<32, 0, true, flip, false>
            <<<1, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
            ((uint4*)keys, (uint4*)keys, numElements, 1 );
    }
    else
    {
        radixSortBlocksKeysOnly<32, 0, false, flip, false>
            <<<1, SORT_CTA_SIZE, 4 * SORT_CTA_SIZE * sizeof(uint)>>>
            ((uint4*)keys, (uint4*)keys, numElements, 1 );
    }

    if (flip)
        unflipFloats<<<1, SORT_CTA_SIZE>>>(keys, numElements);


    CUT_CHECK_ERROR("radixSortSingleBlock");
}

/** 
 * @brief Main radix sort function. For keys only configuration.
 *
 * Main radix sort function.  Sorts in place in the keys array,
 * but uses the other device arrays as temporary storage.  All pointer 
 * parameters are device pointers.  Uses scan for the prefix sum of
 * radix counters.
 * 
 * @param[in,out] keys Keys to be sorted.
 * @param[in] plan Configuration information for RadixSort.
 * @param[in] flipBits Is set true if key datatype is a float (neg. numbers) 
 *        for special float sorting operations.
 * @param[in] numElements Number of elements in the sort.
 * @param[in] keyBits Number of interesting bits in the key
**/
extern "C"
void radixSortKeysOnly(uint *keys,
                       const CUDPPRadixSortPlan *plan, 
                       bool flipBits, 
                       size_t numElements,
                       int keyBits)
{

    if(numElements <= WARP_SIZE)
    {
        if (flipBits)
            radixSortSingleWarpKeysOnly<true><<<1, numElements>>>(keys, numElements);
        else
            radixSortSingleWarpKeysOnly<false><<<1, numElements>>>(keys, numElements);
        return;
    }
    if(numElements <= SORT_CTA_SIZE * 4)
    {
        if (flipBits)
            radixSortSingleBlockKeysOnly<true>(keys, numElements);
        else
            radixSortSingleBlockKeysOnly<false>(keys, numElements);
        return;
    }

    // flip float bits on the first pass, unflip on the last pass
    if (flipBits) 
    {
        radixSortStepKeysOnly<4,  0, true, false>(keys, plan, numElements);
    }
    else
    {
        radixSortStepKeysOnly<4,  0, false, false>(keys, plan, numElements);
    }

    if (keyBits > 4)
    {
        radixSortStepKeysOnly<4,  4, false, false>(keys, plan, numElements);
    }
    if (keyBits > 8)
    {
        radixSortStepKeysOnly<4,  8, false, false>(keys, plan, numElements);
    }
    if (keyBits > 12)
    {
        radixSortStepKeysOnly<4, 12, false, false>(keys, plan, numElements);
    }
    if (keyBits > 16)
    {
        radixSortStepKeysOnly<4, 16, false, false>(keys, plan, numElements);
    }
    if (keyBits > 20)
    {
        radixSortStepKeysOnly<4, 20, false, false>(keys, plan, numElements);
    }
    if (keyBits > 24)
    {
       radixSortStepKeysOnly<4, 24, false, false>(keys, plan, numElements);
    }
    if (keyBits > 28)
    {
        if (flipBits) // last pass
        {
            radixSortStepKeysOnly<4, 28, false, true>(keys, plan, numElements);
        }
        else
        {
            radixSortStepKeysOnly<4, 28, false, false>(keys, plan, numElements);
        }
    }
}

/**
 * @brief Wrapper to call main radix sort function. For floats and keys only.
 *
 * Calls the radixSortKeysOnly function setting parameters for floats.
 * 
 * @param[in,out] keys Keys to be sorted.
 * @param[in] plan Configuration information for RadixSort.
 * @param[in] negativeKeys Is set true if key flipBits is to be true in 
 *                     radixSortKeysOnly().
 * @param[in] numElements Number of elements in the sort.
 * @param[in] keyBits Number of interesting bits in the key
**/
extern "C"
void radixSortFloatKeysOnly(float *keys, 
                            const CUDPPRadixSortPlan *plan,                        
                            bool  negativeKeys,
                            size_t numElements,
                            int keyBits)
{
    radixSortKeysOnly((uint*)keys, plan, negativeKeys, numElements, keyBits);
}

extern "C"
void initDeviceParameters(CUDPPRadixSortPlan *plan)
{
    int deviceID = -1;
    if (cudaSuccess == cudaGetDevice(&deviceID))
    {
        cudaDeviceProp devprop;
        cudaGetDeviceProperties(&devprop, deviceID);

        int smVersion = devprop.major * 10 + devprop.minor;

        // sm_12 and later devices don't need help with coalesce in reorderData kernel
        plan->m_bManualCoalesce = (smVersion < 12);

        // sm_20 and later devices are better off not using persistent CTAs
        plan->m_bUsePersistentCTAs = (smVersion < 20);

        if (plan->m_bUsePersistentCTAs)
        {
            // The following is only true on pre-sm_20 devices (pre-Fermi):
            // Empirically we have found that for some (usually larger) sort
            // sizes it is better to use exactly as many "persistent" CTAs 
            // as can fill the GPU, which loop over the "blocks" of work. For smaller 
            // arrays it is better to use the typical CUDA approach of launching one CTA
            // per block of work.
            // 0-element of these two-element arrays is for key-value sorts
            // 1-element is for key-only sorts
            plan->m_persistentCTAThreshold[0] = plan->m_bManualCoalesce ? 16777216 : 524288;
            plan->m_persistentCTAThresholdFullBlocks[0] = plan->m_bManualCoalesce ? 2097152: 524288;
            plan->m_persistentCTAThreshold[1] = plan->m_bManualCoalesce ? 16777216 : 8388608;
            plan->m_persistentCTAThresholdFullBlocks[1] = plan->m_bManualCoalesce ? 2097152: 0;

            // create a map of function pointers to register counts for more accurate occupancy calculation
            // Must pass in the dynamic shared memory used by each kernel, since the runtime doesn't know it
            // Note we only insert the "loop" version of the kernels (the one with the last template param = true)
            // Because those are the only ones that require persistent CTAs that maximally fill the device.
            computeNumCTAs(radixSortBlocks<4, 0, false, false, true>,         4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            computeNumCTAs(radixSortBlocks<4, 0, false, true,  true>,         4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            computeNumCTAs(radixSortBlocks<4, 0, true, false,  true>,         4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            computeNumCTAs(radixSortBlocks<4, 0, true, true,  true>,          4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            
            computeNumCTAs(radixSortBlocksKeysOnly<4, 0, false, false, true>, 4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            computeNumCTAs(radixSortBlocksKeysOnly<4, 0, false, true, true>,  4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            computeNumCTAs(radixSortBlocksKeysOnly<4, 0, true, false, true>,  4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            computeNumCTAs(radixSortBlocksKeysOnly<4, 0, true, true, true>,   4 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);

            computeNumCTAs(findRadixOffsets<0, false, true>,                  3 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);
            computeNumCTAs(findRadixOffsets<0, true, true>,                   3 * SORT_CTA_SIZE * sizeof(uint), SORT_CTA_SIZE);

            computeNumCTAs(reorderData<0, false, false, false, true>,         0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderData<0, false, false, true, true>,          0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderData<0, false, true, false, true>,          0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderData<0, false, true, true, true>,           0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderData<0, true, false, false, true>,          0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderData<0, true, false, true, true>,           0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderData<0, true, true, false, true>,           0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderData<0, true, true, true, true>,            0,                                SORT_CTA_SIZE);

            computeNumCTAs(reorderDataKeysOnly<0, false, false, false, true>, 0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderDataKeysOnly<0, false, false, true, true>,  0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderDataKeysOnly<0, false, true, false, true>,  0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderDataKeysOnly<0, false, true, true, true>,   0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderDataKeysOnly<0, true, false, false, true>,  0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderDataKeysOnly<0, true, false, true, true>,   0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderDataKeysOnly<0, true, true, false, true>,   0,                                SORT_CTA_SIZE);
            computeNumCTAs(reorderDataKeysOnly<0, true, true, true, true>,    0,                                SORT_CTA_SIZE);
                   
            computeNumCTAs(emptyKernel,                                       0,                                SORT_CTA_SIZE);
        }
    }
}

/**
 * @brief From the programmer-specified sort configuration, 
 *        creates internal memory for performing the sort.
 * 
 * @param[in] plan Pointer to CUDPPRadixSortPlan object
**/
extern "C"
void allocRadixSortStorage(CUDPPRadixSortPlan *plan)
{               
        
    unsigned int numElements = plan->m_numElements;

    unsigned int numBlocks = 
        ((numElements % (SORT_CTA_SIZE * 4)) == 0) ? 
            (numElements / (SORT_CTA_SIZE * 4)) : 
            (numElements / (SORT_CTA_SIZE * 4) + 1);
                        
    switch(plan->m_config.datatype)
    {
    case CUDPP_UINT:
        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_tempKeys, 
                                  numElements * sizeof(unsigned int)));

        if (!plan->m_bKeysOnly)
            CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_tempValues, 
                           numElements * sizeof(unsigned int)));

        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_counters, 
                       WARP_SIZE * numBlocks * sizeof(unsigned int)));

        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_countersSum,
                       WARP_SIZE * numBlocks * sizeof(unsigned int)));

        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_blockOffsets, 
                       WARP_SIZE * numBlocks * sizeof(unsigned int)));
    break;

    case CUDPP_FLOAT:
        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_tempKeys,
                                   numElements * sizeof(float)));

        if (!plan->m_bKeysOnly)
            CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_tempValues,
                           numElements * sizeof(float)));

        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_counters,
                       WARP_SIZE * numBlocks * sizeof(float)));

        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_countersSum,
                       WARP_SIZE * numBlocks * sizeof(float)));

        CUDA_SAFE_CALL(cudaMalloc((void **)&plan->m_blockOffsets,
                       WARP_SIZE * numBlocks * sizeof(float)));     
    break;
    }
        
    initDeviceParameters(plan);
}

/** @brief Deallocates intermediate memory from allocRadixSortStorage.
 *
 *
 * @param[in] plan Pointer to CUDPPRadixSortPlan object
**/
extern "C"
void freeRadixSortStorage(CUDPPRadixSortPlan* plan)
{
    CUDA_SAFE_CALL( cudaFree(plan->m_tempKeys));
    CUDA_SAFE_CALL( cudaFree(plan->m_tempValues));
    CUDA_SAFE_CALL( cudaFree(plan->m_counters));
    CUDA_SAFE_CALL( cudaFree(plan->m_countersSum));
    CUDA_SAFE_CALL( cudaFree(plan->m_blockOffsets));
}

/** @brief Dispatch function to perform a sort on an array with 
 * a specified configuration.
 *
 * This is the dispatch routine which calls radixSort...() with 
 * appropriate template parameters and arguments as specified by 
 * the plan.
 * @param[in,out] keys Keys to be sorted.
 * @param[in,out] values Associated values to be sorted (through keys).
 * @param[in] numElements Number of elements in the sort.
 * @param[in] keyBits Number of interesting bits in the key*
 * @param[in] plan Configuration information for RadixSort.
**/
extern "C"
void cudppRadixSortDispatch(void  *keys,
                            void  *values,
                            size_t numElements,
                            int   keyBits,
                            const CUDPPRadixSortPlan *plan)
{              
    if(plan->m_bKeysOnly)
    {
        switch(plan->m_config.datatype)
        {
        case CUDPP_UINT:
            radixSortKeysOnly((uint*)keys, plan, false, 
                              numElements, keyBits);
            break;
        case CUDPP_FLOAT:
            radixSortFloatKeysOnly((float*)keys, plan, true,
                                    numElements, keyBits);
        }
    }
    else
    {
        switch(plan->m_config.datatype)
        {
        case CUDPP_UINT:      
            radixSort((uint*)keys, (uint*) values, plan, 
                      numElements, false, keyBits);
            break;
        case CUDPP_FLOAT: 
            radixSortFloatKeys((float*)keys, (uint*) values, plan, 
                               numElements, true, keyBits);
        }
    }
}                            

/** @} */ // end radixsort functions
/** @} */ // end cudpp_app
