// -------------------------------------------------------------
// CUDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision$
// $Date$
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt 
// in the root directory of this source distribution.
// ------------------------------------------------------------- 

#include "cudpp_radixsort.h"
#include <cudpp_globals.h>
#include "sharedmem.h"
#include "cta/radixsort_cta.cu"

#ifdef __DEVICE_EMULATION__
#define __EMUSYNC  __syncthreads()
#else
#define __EMUSYNC
#endif

/**
 * @file
 * radixsort_app.cu
 *   
 * @brief CUDPP kernel-level radix sorting routines
 */

/** \addtogroup cudpp_kernel
  * @{
 */

/** @name RadixSort Functions
 * @{
 */



typedef unsigned int uint;

/** @brief And empty kernel used to reset CTA issue hardware
 **/
__global__ void emptyKernel() {}


/** @brief Does special binary arithmetic before sorting floats
 * 
 * Uses floatFlip function to flip bits.
 * @param[in,out] values  Values to be manipulated
 * @param[in] numValues Number of values to be flipped 
 **/

__global__ void 
LAUNCH_BOUNDS(SORT_CTA_SIZE)
flipFloats(uint *values, uint numValues)
{
    uint index = __umul24(blockDim.x*4, blockIdx.x) + threadIdx.x; 
    if (index < numValues) values[index] = floatFlip<true>(values[index]);
    index += blockDim.x;
    if (index < numValues) values[index] = floatFlip<true>(values[index]);
    index += blockDim.x;
    if (index < numValues) values[index] = floatFlip<true>(values[index]);
    index += blockDim.x;
    if (index < numValues) values[index] = floatFlip<true>(values[index]);
}

/** @brief Undoes the flips from flipFloats
 * 
 * Uses floatUnflip function to unflip bits.
 * @param[in,out] values  Values to be manipulated
 * @param[in] numValues Number of values to be unflipped 
 **/
__global__ void 
LAUNCH_BOUNDS(SORT_CTA_SIZE)
unflipFloats(uint *values, uint numValues)
{
    uint index = __umul24(blockDim.x*4, blockIdx.x) + threadIdx.x; 
    if (index < numValues) values[index] = floatUnflip<true>(values[index]);
    index += blockDim.x;
    if (index < numValues) values[index] = floatUnflip<true>(values[index]);
    index += blockDim.x;
    if (index < numValues) values[index] = floatUnflip<true>(values[index]);
    index += blockDim.x;
    if (index < numValues) values[index] = floatUnflip<true>(values[index]);
}


/** @brief Optimization for sorts of WARP_SIZE or fewer elements
 * 
 * @param[in,out] keys  Keys to be sorted.
 * @param[in,out] values Associated values to be sorted (through keys).
 * @param[in] numElements Number of elements in the sort.
 */
template <bool flip>
__global__ 
LAUNCH_BOUNDS(WARP_SIZE)
void radixSortSingleWarp(uint *keys, 
                         uint *values, 
                         uint numElements)
{
    volatile __shared__ uint sKeys[WARP_SIZE]; //remove class distinctions
    volatile __shared__ uint sValues[WARP_SIZE];
    volatile __shared__ uint sFlags[WARP_SIZE];

    sKeys[threadIdx.x]   = floatFlip<flip>(keys[threadIdx.x]);
    sValues[threadIdx.x] = values[threadIdx.x];
    
    __EMUSYNC; // emulation only

    for(uint i = 1; i < numElements; i++)
    {
        uint key_i = sKeys[i];
        uint val_i = sValues[i];
        
        sFlags[threadIdx.x] = 0;
      
        uint temp, tempval;
        if( (threadIdx.x < i) && (sKeys[threadIdx.x] > key_i) ) 
        {
            temp = sKeys[threadIdx.x];
            tempval = sValues[threadIdx.x];
            sFlags[threadIdx.x] = 1;

#ifdef __DEVICE_EMULATION__
        }
        __EMUSYNC;
        if( (threadIdx.x < i) && (sKeys[threadIdx.x] > key_i) ) 
        {
#endif
            sKeys[threadIdx.x + 1] = temp;
            sValues[threadIdx.x + 1] = tempval;
            sFlags[threadIdx.x + 1] = 0;
        }

        
        if(sFlags[threadIdx.x] == 1 )
        {
            sKeys[threadIdx.x] = key_i;
            sValues[threadIdx.x] = val_i;
        }

        __EMUSYNC; // emulation only

    }
    keys[threadIdx.x]   = floatUnflip<flip>(sKeys[threadIdx.x]);
    values[threadIdx.x] = sValues[threadIdx.x];
}


/** @brief Optimization for sorts of WARP_SIZE or fewer elements. Keys-Only version.
 *
 * @param[in,out] keys Keys to be sorted
 * @param[in] numElements Total number of elements to be sorted
**/

template <bool flip>
__global__ 
LAUNCH_BOUNDS(WARP_SIZE)
void radixSortSingleWarpKeysOnly(uint *keys, 
                                 uint numElements)
{
    volatile __shared__ uint sKeys[WARP_SIZE];
    volatile __shared__ uint sFlags[WARP_SIZE];

    sKeys[threadIdx.x]   = floatFlip<flip>(keys[threadIdx.x]);
    
    __EMUSYNC; // emulation only

    for(uint i = 1; i < numElements; i++)
    {
        uint key_i = sKeys[i];
        
        sFlags[threadIdx.x] = 0;
        
        uint temp;
        if( (threadIdx.x < i) && (sKeys[threadIdx.x] > key_i) ) 
        {
            temp = sKeys[threadIdx.x];
            sFlags[threadIdx.x] = 1;
#ifdef __DEVICE_EMULATION__
        }
        __EMUSYNC;
        if( (threadIdx.x < i) && (sKeys[threadIdx.x] > key_i) ) 
        {
#endif
            sKeys[threadIdx.x + 1] = temp;
            sFlags[threadIdx.x + 1] = 0;
        }
        if(sFlags[threadIdx.x] == 1 )
        {
            sKeys[threadIdx.x] = key_i;
        }

        __EMUSYNC; // emulation only

    }
    keys[threadIdx.x]   = floatUnflip<flip>(sKeys[threadIdx.x]);
}

/** @brief sorts all blocks of data independently in shared memory.  
* Each thread block (CTA) sorts one block of 4*CTA_SIZE elements
* 
* The radix sort is done in two stages.  This stage calls radixSortBlock on each 
* block independently, sorting on the basis of bits (startbit) -> (startbit + nbits)
* 
* Template parameters are used to generate efficient code for various special cases
* For example, we have to handle arrays that are a multiple of the block size (fullBlocks)
* differently than arrays that are not.  "flip" is used to only compile in the
* float flip code when float keys are used.  "loop" is used when persistent CTAs
* are used. 
*
* By persistent CTAs we mean that we launch only as many thread blocks as can 
* be resident in the GPU and no more, rather than launching as many threads as
* we have elements. Persistent CTAs loop over blocks of elements until all work
* is complete.  This can be faster in some cases.  In our tests it is faster
* for large sorts (and the threshold is higher on compute version 1.1 and earlier
* GPUs than it is on compute version 1.2 GPUs.
* 
* @param[out] keysOut Output of sorted keys 
* @param[out] valuesOut Output of associated values 
* @param[in]  keysIn Input of unsorted keys in GPU 
* @param[in]  valuesIn Input of associated input values 
* @param[in]  numElements Total number of elements to sort
* @param[in]  totalBlocks The number of blocks of data to sort
*/
template<uint nbits, uint startbit, bool fullBlocks, bool flip, bool loop>
__global__ void 
LAUNCH_BOUNDS(SORT_CTA_SIZE)
radixSortBlocks(uint4* keysOut, uint4* valuesOut, 
                                uint4* keysIn, uint4* valuesIn, 
                                uint numElements, uint totalBlocks)
{
    extern __shared__ uint4 sMem[];

    uint4 key, value;


    uint blockId = blockIdx.x;

    while (!loop || blockId < totalBlocks)
    {
        uint i = blockId * blockDim.x + threadIdx.x;
        uint idx = i << 2;

        // handle non-full last block if array is not multiple of 1024 numElements
        if (!fullBlocks && idx+3 >= numElements)
        {
            if (idx >= numElements)
            {
                key   = make_uint4(UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX);
                value = make_uint4(UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX);
            }
            else
            {
                // for non-full block, we handle uint1 values instead of uint4
                uint *keys1    = (uint*)keysIn;
                uint *values1  = (uint*)valuesIn;

                key.x = (idx   < numElements) ? floatFlip<flip>(keys1[idx])   : UINT_MAX;
                key.y = (idx+1 < numElements) ? floatFlip<flip>(keys1[idx+1]) : UINT_MAX;
                key.z = (idx+2 < numElements) ? floatFlip<flip>(keys1[idx+2]) : UINT_MAX;
                key.w = UINT_MAX;

                value.x = (idx   < numElements) ? values1[idx]   : UINT_MAX;
                value.y = (idx+1 < numElements) ? values1[idx+1] : UINT_MAX;
                value.z = (idx+2 < numElements) ? values1[idx+2] : UINT_MAX;
                value.w = UINT_MAX;
            }
        }
        else
        {
            key = keysIn[i];
            value = valuesIn[i];

            if (flip)
            {
                key.x = floatFlip<flip>(key.x);
                key.y = floatFlip<flip>(key.y);
                key.z = floatFlip<flip>(key.z);
                key.w = floatFlip<flip>(key.w);
            }
        }
        __syncthreads();
        radixSortBlock<nbits, startbit>(key, value);

        // handle non-full last block if array is not multiple of 1024 numElements
        if(!fullBlocks && idx+3 >= numElements)
        {
            if (idx < numElements) 
            {
                // for non-full block, we handle uint1 values instead of uint4
                uint *keys1   = (uint*)keysOut;
                uint *values1 = (uint*)valuesOut;

                keys1[idx]   = key.x;
                values1[idx] = value.x;

                if (idx + 1 < numElements)
                {
                    keys1[idx + 1]   = key.y;
                    values1[idx + 1] = value.y;

                    if (idx + 2 < numElements)
                    {
                        keys1[idx + 2]   = key.z;
                        values1[idx + 2] = value.z;
                    }
                }
            }
        }
        else
        {
            keysOut[i]   = key;
            valuesOut[i] = value;
        }

        if (loop)        
            blockId += gridDim.x;
        else
            break;            
    }
}

/** @brief Computes the number of keys of each radix in each block stores offset.
*
* Given an array with blocks sorted according to a 4-bit radix group, each 
* block counts the number of keys that fall into each radix in the group, and 
* finds the starting offset of each radix in the block.  It then writes the radix 
* counts to the counters array, and the starting offsets to the blockOffsets array.
*
* Template parameters are used to generate efficient code for various special cases
* For example, we have to handle arrays that are a multiple of the block size 
* (fullBlocks) differently than arrays that are not. "loop" is used when persistent 
* CTAs are used. 
*
* By persistent CTAs we mean that we launch only as many thread blocks as can 
* be resident in the GPU and no more, rather than launching as many threads as
* we have elements. Persistent CTAs loop over blocks of elements until all work
* is complete.  This can be faster in some cases.  In our tests it is faster
* for large sorts (and the threshold is higher on compute version 1.1 and earlier
* GPUs than it is on compute version 1.2 GPUs.
* 
* @param[in] keys Input keys
* @param[out] counters Radix count for each block
* @param[out] blockOffsets The offset address for each block
* @param[in] numElements Total number of elements
* @param[in] totalBlocks Total number of blocks
**/
template<uint startbit, bool fullBlocks, bool loop>
__global__ void 
LAUNCH_BOUNDS(SORT_CTA_SIZE)
findRadixOffsets(uint2 *keys, 
                 uint  *counters, 
                 uint  *blockOffsets, 
                 uint   numElements,
                 uint   totalBlocks)
{
    extern __shared__ uint sRadix1[];
    __shared__ uint  sStartPointers[16];

    uint blockId = blockIdx.x;   

    while (!loop || blockId < totalBlocks)
    {
        uint2 radix2;

        uint i       = blockId * blockDim.x + threadIdx.x;

        // handle non-full last block if array is not multiple of 1024 numElements
        if(!fullBlocks && ((i + 1) << 1 ) > numElements )
        {
            // handle uint1 rather than uint2 for non-full blocks
            uint *keys1 = (uint*)keys;
            uint j = i << 1; 

            radix2.x = (j < numElements) ? keys1[j] : UINT_MAX; 
            j++;
            radix2.y = (j < numElements) ? keys1[j] : UINT_MAX;
        }
        else
        {
            radix2 = keys[i];
        }

        sRadix1[2 * threadIdx.x]     = (radix2.x >> startbit) & 0xF;
        sRadix1[2 * threadIdx.x + 1] = (radix2.y >> startbit) & 0xF;

        // Finds the position where the sRadix1 entries differ and stores start 
        // index for each radix.
        if(threadIdx.x < 16) 
        { 
            sStartPointers[threadIdx.x] = 0; 
        }
        __syncthreads();

        if((threadIdx.x > 0) && (sRadix1[threadIdx.x] != sRadix1[threadIdx.x - 1]) ) 
        {
            sStartPointers[sRadix1[threadIdx.x]] = threadIdx.x;
        }
        if(sRadix1[threadIdx.x + SORT_CTA_SIZE] != sRadix1[threadIdx.x + SORT_CTA_SIZE - 1]) 
        {
            sStartPointers[sRadix1[threadIdx.x + SORT_CTA_SIZE]] = threadIdx.x + SORT_CTA_SIZE;
        }
        __syncthreads();

        if(threadIdx.x < 16) 
        {
            blockOffsets[blockId*16 + threadIdx.x] = sStartPointers[threadIdx.x];
        }
        __syncthreads();

        // Compute the sizes of each block.
        if((threadIdx.x > 0) && (sRadix1[threadIdx.x] != sRadix1[threadIdx.x - 1]) ) 
        {
            sStartPointers[sRadix1[threadIdx.x - 1]] = 
                threadIdx.x - sStartPointers[sRadix1[threadIdx.x - 1]];
        }
        if(sRadix1[threadIdx.x + SORT_CTA_SIZE] != sRadix1[threadIdx.x + SORT_CTA_SIZE - 1] ) 
        {
            sStartPointers[sRadix1[threadIdx.x + SORT_CTA_SIZE - 1]] = 
                threadIdx.x + SORT_CTA_SIZE - sStartPointers[sRadix1[threadIdx.x + SORT_CTA_SIZE - 1]];
        }


        if(threadIdx.x == SORT_CTA_SIZE - 1) 
        {
            sStartPointers[sRadix1[2 * SORT_CTA_SIZE - 1]] = 
                2 * SORT_CTA_SIZE - sStartPointers[sRadix1[2 * SORT_CTA_SIZE - 1]];
        }
        __syncthreads();

        if(threadIdx.x < 16) 
        {
            counters[threadIdx.x * totalBlocks + blockId] = 
                sStartPointers[threadIdx.x];
        }

        if (loop)
            blockId += gridDim.x;
        else
            break;
    }
}


/**@brief Reorders data in the global array.
*
* reorderData shuffles data in the array globally after the radix
* offsets have been found. On compute version 1.1 and earlier GPUs, this code depends 
* on SORT_CTA_SIZE being 16 * number of radices (i.e. 16 * 2^nbits).
* 
* On compute version 1.1 GPUs ("manualCoalesce=true") this function ensures
* that all writes are coalesced using extra work in the kernel.  On later
* GPUs coalescing rules have been relaxed, so this extra overhead hurts 
* performance.  On these GPUs we set manualCoalesce=false and directly store
* the results.
*
* Template parameters are used to generate efficient code for various special cases
* For example, we have to handle arrays that are a multiple of the block size 
* (fullBlocks) differently than arrays that are not.  "loop" is used when persistent 
* CTAs are used. 
*
* By persistent CTAs we mean that we launch only as many thread blocks as can 
* be resident in the GPU and no more, rather than launching as many threads as
* we have elements. Persistent CTAs loop over blocks of elements until all work
* is complete.  This can be faster in some cases.  In our tests it is faster
* for large sorts (and the threshold is higher on compute version 1.1 and earlier
* GPUs than it is on compute version 1.2 GPUs.
*
* @param[out] outKeys Output of sorted keys 
* @param[out] outValues Output of associated values 
* @param[in] keys Input of unsorted keys in GPU 
* @param[in] values Input of associated input values 
* @param[in] blockOffsets The offset address for each block
* @param[in] offsets Address of each radix within each block
* @param[in] sizes Number of elements in a block
* @param[in] numElements Total number of elements
* @param[in] totalBlocks Total number of data blocks to process
*
* @todo Args that are const below should be prototyped as const
**/
template<uint startbit, bool fullBlocks, bool manualCoalesce, bool unflip, bool loop>
__global__ void 
LAUNCH_BOUNDS(SORT_CTA_SIZE)
reorderData(uint  *outKeys, 
            uint  *outValues, 
            uint2 *keys, 
            uint2 *values, 
            uint  *blockOffsets, 
            uint  *offsets, 
            uint  *sizes, 
            uint   numElements,
            uint   totalBlocks)
{
    __shared__ uint2 sKeys2[SORT_CTA_SIZE];
    __shared__ uint2 sValues2[SORT_CTA_SIZE];
    __shared__ uint sOffsets[16];
    __shared__ uint sBlockOffsets[16];

    uint *sKeys1   = (uint*)sKeys2; 
    uint *sValues1 = (uint*)sValues2; 

    uint blockId = blockIdx.x;   

    while (!loop || blockId < totalBlocks)
    {
        uint i = blockId * blockDim.x + threadIdx.x;

        // handle non-full last block if array is not multiple of 1024 numElements
        if(!fullBlocks && (((i + 1) << 1) > numElements))
        {
            uint *keys1   = (uint*)keys;
            uint *values1 = (uint*)values;
            uint j = i << 1; 

            sKeys1[threadIdx.x << 1]   = (j < numElements) ? keys1[j]   : UINT_MAX; 
            sValues1[threadIdx.x << 1] = (j < numElements) ? values1[j] : UINT_MAX; 
            j++; 
            sKeys1[(threadIdx.x << 1) + 1]   = (j < numElements) ? keys1[j]   : UINT_MAX; 
            sValues1[(threadIdx.x << 1) + 1] = (j < numElements) ? values1[j] : UINT_MAX; 
        }
        else
        {
            sKeys2[threadIdx.x]   = keys[i];
            sValues2[threadIdx.x] = values[i];
        }

        if (!manualCoalesce)
        {
            if(threadIdx.x < 16)  
            {
                sOffsets[threadIdx.x]      = offsets[threadIdx.x * totalBlocks + blockId];
                sBlockOffsets[threadIdx.x] = blockOffsets[blockId * 16 + threadIdx.x];
            }
            __syncthreads();

            uint radix = (sKeys1[threadIdx.x] >> startbit) & 0xF;
            uint globalOffset = sOffsets[radix] + threadIdx.x - sBlockOffsets[radix];

            if (fullBlocks || globalOffset < numElements)
            {
                outKeys[globalOffset]   = floatUnflip<unflip>(sKeys1[threadIdx.x]);
                outValues[globalOffset] = sValues1[threadIdx.x];
            }

            radix = (sKeys1[threadIdx.x + SORT_CTA_SIZE] >> startbit) & 0xF;
            globalOffset = sOffsets[radix] + threadIdx.x + SORT_CTA_SIZE - sBlockOffsets[radix];

            if (fullBlocks || globalOffset < numElements)
            {
                outKeys[globalOffset]   = floatUnflip<unflip>(sKeys1[threadIdx.x + SORT_CTA_SIZE]);
                outValues[globalOffset] = sValues1[threadIdx.x + SORT_CTA_SIZE];
            }
        }
        else
        {
            __shared__ uint sSizes[16];

            if(threadIdx.x < 16)  
            {
                sOffsets[threadIdx.x]      = offsets[threadIdx.x * totalBlocks + blockId];
                sBlockOffsets[threadIdx.x] = blockOffsets[blockId * 16 + threadIdx.x];
                sSizes[threadIdx.x]        = sizes[threadIdx.x * totalBlocks + blockId];
            }
            __syncthreads();

            // 1 half-warp is responsible for writing out all values for 1 radix. 
            // Loops if there are more than 16 values to be written out. 
            // All start indices are rounded down to the nearest multiple of 16, and
            // all end indices are rounded up to the nearest multiple of 16.
            // Thus it can do extra work if the start and end indices are not multiples of 16
            // This is bounded by a factor of 2 (it can do 2X more work at most).

            const uint halfWarpID     = threadIdx.x >> 4;

            const uint halfWarpOffset = threadIdx.x & 0xF;
            const uint leadingInvalid = sOffsets[halfWarpID] & 0xF;

            uint startPos = sOffsets[halfWarpID] & 0xFFFFFFF0;
            uint endPos   = (sOffsets[halfWarpID] + sSizes[halfWarpID]) + 15 - 
                ((sOffsets[halfWarpID] + sSizes[halfWarpID] - 1) & 0xF);
            uint numIterations = endPos - startPos;

            uint outOffset = startPos + halfWarpOffset;
            uint inOffset  = sBlockOffsets[halfWarpID] - leadingInvalid + halfWarpOffset;

            for(uint j = 0; j < numIterations; j += 16, outOffset += 16, inOffset += 16)
            {       
                if( (outOffset >= sOffsets[halfWarpID]) && 
                    (inOffset - sBlockOffsets[halfWarpID] < sSizes[halfWarpID])) 
                {
                    if(blockId < totalBlocks - 1 || outOffset < numElements) 
                    {
                        outKeys[outOffset]   = floatUnflip<unflip>(sKeys1[inOffset]);
                        outValues[outOffset] = sValues1[inOffset];
                    }
                }       
            }
        }

        if (loop)
        {
            blockId += gridDim.x;
            __syncthreads();
        }
        else
            break;
    }
}

/** @brief Sorts all blocks of data independently in shared memory.  
*  Each thread block (CTA) sorts one block of 4*CTA_SIZE elements
* 
* The radix sort is done in two stages.  This stage calls radixSortBlock on each 
* block independently, sorting on the basis of bits (startbit) -> (startbit + nbits)
* 
* Template parameters are used to generate efficient code for various special cases
* For example, we have to handle arrays that are a multiple of the block size (fullBlocks)
* differently than arrays that are not.  "flip" is used to only compile in the
* float flip code when float keys are used.  "loop" is used when persistent CTAs
* are used. 
*
* By persistent CTAs we mean that we launch only as many thread blocks as can 
* be resident in the GPU and no more, rather than launching as many threads as
* we have elements. Persistent CTAs loop over blocks of elements until all work
* is complete.  This can be faster in some cases.  In our tests it is faster
* for large sorts (and the threshold is higher on compute version 1.1 and earlier
* GPUs than it is on compute version 1.2 GPUs.
* 
* @param[out] keysOut Output of sorted keys GPU main memory
* @param[in] keysIn Input of unsorted keys in GPU main memory
* @param[in] numElements Total number of elements to sort
* @param[in] totalBlocks Total number of blocks to sort
*
*/
template<uint nbits, uint startbit, bool fullBlocks, bool flip, bool loop>
__global__ void 
LAUNCH_BOUNDS(SORT_CTA_SIZE)
radixSortBlocksKeysOnly(uint4* keysOut, uint4* keysIn, uint numElements, uint totalBlocks)
{
    extern __shared__ uint4 sMem[];

    uint4 key;

    uint blockId = blockIdx.x;

    while (!loop || blockId < totalBlocks)
    {
        uint i = blockId * blockDim.x + threadIdx.x;
        uint idx = i << 2;

        // handle non-full last block if array is not multiple of 1024 numElements
        if (!fullBlocks && idx+3 >= numElements)
        {
            if (idx >= numElements)
            {
                key   = make_uint4(UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX);
            }
            else
            {
                // for non-full block, we handle uint1 values instead of uint4
                uint *keys1    = (uint*)keysIn;

                key.x = (idx   < numElements) ? floatFlip<flip>(keys1[idx])   : UINT_MAX;
                key.y = (idx+1 < numElements) ? floatFlip<flip>(keys1[idx+1]) : UINT_MAX;
                key.z = (idx+2 < numElements) ? floatFlip<flip>(keys1[idx+2]) : UINT_MAX;
                key.w = UINT_MAX;
            }
        }
        else
        {
            key = keysIn[i];
            if (flip)
            {
                key.x = floatFlip<flip>(key.x);
                key.y = floatFlip<flip>(key.y);
                key.z = floatFlip<flip>(key.z);
                key.w = floatFlip<flip>(key.w);
            }            
        }
        __syncthreads();
        radixSortBlockKeysOnly<nbits, startbit>(key);

        // handle non-full last block if array is not multiple of 1024 numElements
        if(!fullBlocks && idx+3 >= numElements)
        {
            if (idx < numElements) 
            {
                // for non-full block, we handle uint1 values instead of uint4
                uint *keys1   = (uint*)keysOut;

                keys1[idx]   = key.x;

                if (idx + 1 < numElements)
                {
                    keys1[idx + 1]   = key.y;

                    if (idx + 2 < numElements)
                    {
                        keys1[idx + 2]   = key.z;
                    }
                }
            }
        }
        else
        {
            keysOut[i]   = key;
        }

        if (loop)
            blockId += gridDim.x;
        else
            break;
    }
}

/** @brief Reorders data in the global array.
*
* reorderDataKeysOnly shuffles data in the array globally after the radix offsets 
* have been found. On compute version 1.1 and earlier GPUs, this code depends 
* on SORT_CTA_SIZE being 16 * number of radices (i.e. 16 * 2^nbits).
* 
* On compute version 1.1 GPUs ("manualCoalesce=true") this function ensures
* that all writes are coalesced using extra work in the kernel.  On later
* GPUs coalescing rules have been relaxed, so this extra overhead hurts 
* performance.  On these GPUs we set manualCoalesce=false and directly store
* the results.
*
* Template parameters are used to generate efficient code for various special cases
* For example, we have to handle arrays that are a multiple of the block size 
* (fullBlocks) differently than arrays that are not.  "loop" is used when persistent 
* CTAs are used. 
*
* By persistent CTAs we mean that we launch only as many thread blocks as can 
* be resident in the GPU and no more, rather than launching as many threads as
* we have elements. Persistent CTAs loop over blocks of elements until all work
* is complete.  This can be faster in some cases.  In our tests it is faster
* for large sorts (and the threshold is higher on compute version 1.1 and earlier
* GPUs than it is on compute version 1.2 GPUs.
* 
* @param[out] outKeys Output result of reorderDataKeysOnly()
* @param[in] keys Keys to be reordered
* @param[in] blockOffsets Start offset for each block
* @param[in] offsets Offset of each radix within each block
* @param[in] sizes Number of elements in a block
* @param[in] numElements Total number of elements
* @param[in] totalBlocks Total number of blocks
*/
template<uint startbit, bool fullBlocks, bool manualCoalesce, bool unflip, bool loop>
__global__ void 
LAUNCH_BOUNDS(SORT_CTA_SIZE)
reorderDataKeysOnly(uint  *outKeys, 
                                    uint2 *keys, 
                                    uint  *blockOffsets, 
                                    uint  *offsets, 
                                    uint  *sizes, 
                                    uint   numElements,
                                    uint   totalBlocks)
{
    __shared__ uint2 sKeys2[SORT_CTA_SIZE];
    __shared__ uint sOffsets[16];
    __shared__ uint sBlockOffsets[16];

    uint *sKeys1   = (uint*)sKeys2; 

    uint blockId = blockIdx.x;

    while (!loop || blockId < totalBlocks)
    {
        uint i = blockId * blockDim.x + threadIdx.x;

        // handle non-full last block if array is not multiple of 1024 numElements
        if(!fullBlocks && (((i + 1) << 1) > numElements))
        {
            uint *keys1   = (uint*)keys;
            uint j = i << 1; 

            sKeys1[threadIdx.x << 1]   = (j < numElements) ? keys1[j]   : UINT_MAX; 
            j++; 
            sKeys1[(threadIdx.x << 1) + 1]   = (j < numElements) ? keys1[j]   : UINT_MAX; 
        }
        else
        {
            sKeys2[threadIdx.x]   = keys[i];
        }

        if (!manualCoalesce)
        {
            if(threadIdx.x < 16)  
            {
                sOffsets[threadIdx.x]      = offsets[threadIdx.x * totalBlocks + blockId];
                sBlockOffsets[threadIdx.x] = blockOffsets[blockId * 16 + threadIdx.x];
            }
            __syncthreads();

            uint radix = (sKeys1[threadIdx.x] >> startbit) & 0xF;
            uint globalOffset = sOffsets[radix] + threadIdx.x - sBlockOffsets[radix];

            if (fullBlocks || globalOffset < numElements)
            {
                outKeys[globalOffset]   = floatUnflip<unflip>(sKeys1[threadIdx.x]);
            }

            radix = (sKeys1[threadIdx.x + SORT_CTA_SIZE] >> startbit) & 0xF;
            globalOffset = sOffsets[radix] + threadIdx.x + SORT_CTA_SIZE - sBlockOffsets[radix];

            if (fullBlocks || globalOffset < numElements)
            {
                outKeys[globalOffset]   = floatUnflip<unflip>(sKeys1[threadIdx.x + SORT_CTA_SIZE]);
            }
        }
        else
        {
            __shared__ uint sSizes[16];

            if(threadIdx.x < 16)  
            {
                sOffsets[threadIdx.x]      = offsets[threadIdx.x * totalBlocks + blockId];
                sBlockOffsets[threadIdx.x] = blockOffsets[blockId * 16 + threadIdx.x];
                sSizes[threadIdx.x]        = sizes[threadIdx.x * totalBlocks + blockId];
            }
            __syncthreads();

            // 1 half-warp is responsible for writing out all values for 1 radix. 
            // Loops if there are more than 16 values to be written out. 
            // All start indices are rounded down to the nearest multiple of 16, and
            // all end indices are rounded up to the nearest multiple of 16.
            // Thus it can do extra work if the start and end indices are not multiples of 16
            // This is bounded by a factor of 2 (it can do 2X more work at most).

            const uint halfWarpID     = threadIdx.x >> 4;

            const uint halfWarpOffset = threadIdx.x & 0xF;
            const uint leadingInvalid = sOffsets[halfWarpID] & 0xF;

            uint startPos = sOffsets[halfWarpID] & 0xFFFFFFF0;
            uint endPos   = (sOffsets[halfWarpID] + sSizes[halfWarpID]) + 15 - 
                ((sOffsets[halfWarpID] + sSizes[halfWarpID] - 1) & 0xF);
            uint numIterations = endPos - startPos;

            uint outOffset = startPos + halfWarpOffset;
            uint inOffset  = sBlockOffsets[halfWarpID] - leadingInvalid + halfWarpOffset;

            for(uint j = 0; j < numIterations; j += 16, outOffset += 16, inOffset += 16)
            {       
                if( (outOffset >= sOffsets[halfWarpID]) && 
                    (inOffset - sBlockOffsets[halfWarpID] < sSizes[halfWarpID])) 
                {
                    if(blockId < totalBlocks - 1 || outOffset < numElements) 
                    {
                        outKeys[outOffset] = floatUnflip<unflip>(sKeys1[inOffset]);
                    }
                }       
            }
        }

        if (loop)
        {
            blockId += gridDim.x;
            __syncthreads();
        }
        else
            break;
    }
}

/** @} */ // end radixsort functions
/** @} */ // end cudpp_kernel
