// -------------------------------------------------------------
// CUDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision$
// $Date$
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt 
// in the root directory of this source distribution.
// ------------------------------------------------------------- 
#include <cudpp_globals.h>
#include "cudpp_radixsort.h"
#include "cta/scan_cta.cu"
#include <cudpp.h>
#include <stdio.h>
 
#include <cudpp_util.h>
#include <math.h>
#include "sharedmem.h"


#ifdef __DEVICE_EMULATION__
#define __EMUSYNC __syncthreads()
#else
#define __EMUSYNC
#endif

/**
 * @file
 * sort_cta.cu
 * 
 * @brief CUDPP CTA-level sort routines
 */

/** \addtogroup cudpp_cta 
* @{
*/

/** @name Radix Sort Functions
* @{
*/


typedef unsigned int uint;

/**
 * @brief Flips bits of single-precision floating-point number (parameterized by doFlip)
 * 
 *  flip a float for sorting
 *  finds SIGN of fp number.
 *  if it's 1 (negative float), it flips all bits
 *  if it's 0 (positive float), it flips the sign only
 * @param[in] f floating-point input (passed as unsigned int)
 * @see floatUnflip
**/

template <bool doFlip>
__device__ uint floatFlip(uint f)
{
    if (doFlip)
    {
        uint mask = -int(f >> 31) | 0x80000000;
        return f ^ mask;
    }
    else
        return f;
}

/**
 * @brief Reverses bit-flip of single-precision floating-point number (parameterized by doFlip)
 * 
 * flip a float back (invert FloatFlip)
 *  signed was flipped from above, so:
 *  if sign is 1 (negative), it flips the sign bit back
 *  if sign is 0 (positive), it flips all bits back
 * @param[in] f floating-point input (passed as unsigned int)
 * @see floatFlip
**/
template <bool doFlip>
__device__ uint floatUnflip(uint f)
{
    if (doFlip)
    {
        uint mask = ((f >> 31) - 1) | 0x80000000;
        return f ^ mask;
    }
    else
        return f;
}

/**
 * @brief Scans one warp quickly, optimized for 32-element warps, using shared memory
 * 
 * Scans each warp in parallel ("warp-scan"), one element per thread.
 * uses 2 numElements of shared memory per thread (64 numElements per warp)
 * 
 * @param[in] val Elements per thread to scan
 * @param[in,out] sData
**/
template<class T, int maxlevel>
__device__ T scanwarp(T val, volatile T* sData)
{
    // The following is the same as 2 * WARP_SIZE * warpId + threadInWarp = 
    // 64*(threadIdx.x >> 5) + (threadIdx.x & (WARP_SIZE - 1))
    int idx = 2 * threadIdx.x - (threadIdx.x & (WARP_SIZE - 1));
    sData[idx] = 0;
    idx += WARP_SIZE;
    T t = sData[idx] = val;          __EMUSYNC;

#ifdef __DEVICE_EMULATION__             
        t = sData[idx -  1]; __EMUSYNC; 
        sData[idx] += t;       __EMUSYNC;
        t = sData[idx -  2];   __EMUSYNC; 
        sData[idx] += t;       __EMUSYNC;
        t = sData[idx -  4];   __EMUSYNC; 
        sData[idx] += t;       __EMUSYNC;
        t = sData[idx -  8];   __EMUSYNC; 
        sData[idx] += t;       __EMUSYNC;
        t = sData[idx - 16];   __EMUSYNC; 
        sData[idx] += t;       __EMUSYNC;
#else
        if (0 <= maxlevel) { sData[idx] = t = t + sData[idx - 1]; } __EMUSYNC;
        if (1 <= maxlevel) { sData[idx] = t = t + sData[idx - 2]; } __EMUSYNC;
        if (2 <= maxlevel) { sData[idx] = t = t + sData[idx - 4]; } __EMUSYNC;
        if (3 <= maxlevel) { sData[idx] = t = t + sData[idx - 8]; } __EMUSYNC;
        if (4 <= maxlevel) { sData[idx] = t = t + sData[idx -16]; } __EMUSYNC;
#endif          
        return sData[idx] - val;  // convert inclusive -> exclusive
}

/**
 * @brief Scans 4*CTA_SIZE unsigned ints in a block
 *
 * scan4 scans 4*CTA_SIZE numElements in a block (4 per
 * thread), using a warp-scan algorithm
 * 
 * @param[in] idata 4-vector of integers to scan
**/
__device__ uint4 scan4(uint4 idata)
{    
    extern  __shared__  uint ptr[];
    
    uint idx = threadIdx.x;

    uint4 val4 = idata;
    uint sum[3];
    sum[0] = val4.x;
    sum[1] = val4.y + sum[0];
    sum[2] = val4.z + sum[1];
    
    uint val = val4.w + sum[2];
    
    val = scanwarp<uint, 4>(val, ptr);
    __syncthreads();

    if ((idx & (WARP_SIZE - 1)) == WARP_SIZE - 1)
    {
        ptr[idx >> 5] = val + val4.w + sum[2];
    }
    __syncthreads();

#ifndef __DEVICE_EMULATION__
    if (idx < WARP_SIZE)
#endif
    {
        ptr[idx] = scanwarp<uint, 2>(ptr[idx], ptr);
    }
    __syncthreads();

    val += ptr[idx >> 5];

    val4.x = val;
    val4.y = val + sum[0];
    val4.z = val + sum[1];
    val4.w = val + sum[2];      
        
    return val4;
}

/**
 * @brief Computes output position for each thread given predicate; trues come first then falses
 * 
 * Rank is the core of the radix sort loop.  Given a predicate, it
 * computes the output position for each thread in an ordering where all
 * True threads come first, followed by all False threads. 
 * This version handles 4 predicates per thread; hence, "rank4".
 *
 * @param[in] preds true/false values for each of the 4 elements in this thread
 *
 * @todo is the description of "preds" correct?
**/
template <int ctasize>
__device__ uint4 rank4(uint4 preds)
{
    uint4 address = scan4(preds);  

    __shared__ uint numtrue;
    if (threadIdx.x == ctasize-1)
    {
        numtrue = address.w + preds.w;
    }
    __syncthreads();

    uint4 rank;
    uint idx = threadIdx.x << 2;
    rank.x = (preds.x) ? address.x : numtrue + idx   - address.x;
    rank.y = (preds.y) ? address.y : numtrue + idx + 1 - address.y;
    rank.z = (preds.z) ? address.z : numtrue + idx + 2 - address.z;
    rank.w = (preds.w) ? address.w : numtrue + idx + 3 - address.w;     
                
    return rank;
}

/**
 * @brief Sorts one block
 *
 * Uses rank to sort one bit at a time: Sorts a block according
 * to bits startbit -> nbits + startbit
 * @param[in,out] key
 * @param[in,out] value
**/
template<uint nbits, uint startbit>
__device__ void radixSortBlock(uint4 &key, uint4 &value)
{
    extern __shared__ uint sMem1[];
    for(uint shift = startbit; shift < (startbit + nbits); ++shift)
    {        
        uint4 lsb;
        lsb.x = !((key.x >> shift) & 0x1);
        lsb.y = !((key.y >> shift) & 0x1);
        lsb.z = !((key.z >> shift) & 0x1);
        lsb.w = !((key.w >> shift) & 0x1); 

        uint4 r = rank4<256>(lsb);

#if 1
        // This arithmetic strides the ranks across 4 SORT_CTA_SIZE regions
        sMem1[(r.x & 3) * SORT_CTA_SIZE + (r.x >> 2)] = key.x;
        sMem1[(r.y & 3) * SORT_CTA_SIZE + (r.y >> 2)] = key.y;
        sMem1[(r.z & 3) * SORT_CTA_SIZE + (r.z >> 2)] = key.z;
        sMem1[(r.w & 3) * SORT_CTA_SIZE + (r.w >> 2)] = key.w; 
        __syncthreads();

        // The above allows us to read without 4-way bank conflicts:
        key.x = sMem1[threadIdx.x];
        key.y = sMem1[threadIdx.x +     SORT_CTA_SIZE];
        key.z = sMem1[threadIdx.x + 2 * SORT_CTA_SIZE];
        key.w = sMem1[threadIdx.x + 3 * SORT_CTA_SIZE];

        __syncthreads();

        sMem1[(r.x & 3) * SORT_CTA_SIZE + (r.x >> 2)] = value.x;
        sMem1[(r.y & 3) * SORT_CTA_SIZE + (r.y >> 2)] = value.y;
        sMem1[(r.z & 3) * SORT_CTA_SIZE + (r.z >> 2)] = value.z;
        sMem1[(r.w & 3) * SORT_CTA_SIZE + (r.w >> 2)] = value.w;
        __syncthreads();

        value.x = sMem1[threadIdx.x];
        value.y = sMem1[threadIdx.x +     SORT_CTA_SIZE];
        value.z = sMem1[threadIdx.x + 2 * SORT_CTA_SIZE];
        value.w = sMem1[threadIdx.x + 3 * SORT_CTA_SIZE];
#else
        sMem1[r.x] = key.x;
        sMem1[r.y] = key.y;
        sMem1[r.z] = key.z;
        sMem1[r.w] = key.w;
        __syncthreads();

        // This access has 4-way bank conflicts
        key = sMem[threadIdx.x];

        __syncthreads();

        sMem1[r.x] = value.x;
        sMem1[r.y] = value.y;
        sMem1[r.z] = value.z;
        sMem1[r.w] = value.w;
        __syncthreads();

        value = sMem[threadIdx.x];
#endif

        __syncthreads();
    }
}

/**
 * @brief Sorts one block. Key-only version.
 *
 * Uses rank to sort one bit at a time: Sorts a block according
 * to bits startbit -> nbits + startbit
 * @param[in,out] key
**/

template<uint nbits, uint startbit>
__device__ void radixSortBlockKeysOnly(uint4 &key)
{
    extern __shared__ uint sMem1[];
    for(uint shift = startbit; shift < (startbit + nbits); ++shift)
    {                   
        uint4 lsb;
        lsb.x = !((key.x >> shift) & 0x1);
        lsb.y = !((key.y >> shift) & 0x1);
        lsb.z = !((key.z >> shift) & 0x1);
        lsb.w = !((key.w >> shift) & 0x1);

        uint4 r = rank4<256>(lsb);

#if 1
        // This arithmetic strides the ranks across 4 CTA_SIZE regions
        sMem1[(r.x & 3) * SORT_CTA_SIZE + (r.x >> 2)] = key.x;
        sMem1[(r.y & 3) * SORT_CTA_SIZE + (r.y >> 2)] = key.y;
        sMem1[(r.z & 3) * SORT_CTA_SIZE + (r.z >> 2)] = key.z;
        sMem1[(r.w & 3) * SORT_CTA_SIZE + (r.w >> 2)] = key.w;
        __syncthreads();

        // The above allows us to read without 4-way bank conflicts:
        key.x = sMem1[threadIdx.x];
        key.y = sMem1[threadIdx.x +     SORT_CTA_SIZE];
        key.z = sMem1[threadIdx.x + 2 * SORT_CTA_SIZE];
        key.w = sMem1[threadIdx.x + 3 * SORT_CTA_SIZE];
#else
        sMem1[r.x] = key.x;
        sMem1[r.y] = key.y;
        sMem1[r.z] = key.z;
        sMem1[r.w] = key.w;
        __syncthreads();

        // This access has 4-way bank conflicts
        key = sMem[threadIdx.x];
#endif

        __syncthreads();
    }
}

/** @} */ // end radix sort functions
/** @} */ // end cudpp_cta
