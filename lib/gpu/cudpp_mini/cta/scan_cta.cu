// ------------------------------------------------------------- 
//  cuDPP -- CUDA Data Parallel Primitives library
//  -------------------------------------------------------------
//  $Revision: 5633 $
//  $Date: 2009-07-01 15:02:51 +1000 (Wed, 01 Jul 2009) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt 
// in the root directory of this source distribution.
// ------------------------------------------------------------- 

/**
 * @file
 * scan_cta.cu
 *
 * @brief CUDPP CTA-level scan routines
 */

/** \defgroup cudpp_cta CUDPP CTA-Level API
  * The CUDPP CTA-Level API contains functions that run on the GPU 
  * device.  These are CUDA \c __device__ functions that are called
  * from within other CUDA device functions (typically 
  * \link cudpp_kernel CUDPP Kernel-Level API\endlink functions).
  * They are called CTA-level functions because they typically process
  * s_data "owned" by each CTA within shared memory, and are agnostic of
  * any other CTAs that may be running (or how many CTAs are running),
  * other than to compute appropriate global memory addresses.
  * @{
  */

/** @name Scan Functions
* @{
*/

#include <cudpp_globals.h>
#include <cudpp_util.h>
#include <math.h>
#include <cudpp.h>

/**
 * @brief Macro to insert necessary __syncthreads() in device emulation mode
 */
#ifdef __DEVICE_EMULATION__
#define __EMUSYNC __syncthreads()
#else
#define __EMUSYNC
#endif

/** 
  * @brief Template class containing compile-time parameters to the scan functions
  *
  * ScanTraits is passed as a template parameter to all scan functions.  By 
  * using these compile-time functions we can enable generic code while 
  * maintaining the highest performance.  This is crucial for the performance 
  * of low-level workhorse algorithms like scan.
  *
  * @param T The datatype of the scan
  * @param oper The ::CUDPPOperator to use for the scan (add, max, etc.)
  * @param multiRow True if this is a multi-row scan
  * @param unroll True if scan inner loops should be unrolled
  * @param sums True if each block should write it's sum to the d_blockSums array (false for single-block scans)
  * @param backward True if this is a backward scan
  * @param fullBlock True if all blocks in this scan are full (CTA_SIZE * SCAN_ELEMENTS_PER_THREAD elements)
  * @param exclusive True for exclusive scans, false for inclusive scans
  */
template <class T, CUDPPOperator oper, bool backward, bool exclusive,
          bool multiRow, bool sums, bool fullBlock>
class ScanTraits
{
public:
    
    //! Returns true if this is a backward scan
    static __device__ bool isBackward()    { return backward; };
    //! Returns true if this is an exclusive scan
    static __device__ bool isExclusive()  { return exclusive; };
    //! Returns true if this a multi-row scan.
    static __device__ bool isMultiRow()    { return multiRow; };
    //! Returns true if this scan writes the sum of each block to the d_blockSums array (multi-block scans)
    static __device__ bool writeSums()     { return sums; };
    //! Returns true if this is a full scan -- all blocks process CTA_SIZE * SCAN_ELEMENTS_PER_THREAD elements
    static __device__ bool isFullBlock()   { return fullBlock; };
    
        
    //! The operator function used for the scan
    static __device__ T op(const T a, const T b)
    {
        return Operator<T, oper>::op(a, b);
    }  

    //! The identity value used by the scan
    static __device__ T identity() { return Operator<T, oper>::identity(); }
};

//! This is used to insert syncthreads to avoid perf loss caused by 128-bit 
//! load overlap that happens on G80.  This gives about a 15% boost on scans on 
//! G80.
//! @todo Parameterize this in case this perf detail changes on future GPUs.
#define DISALLOW_LOADSTORE_OVERLAP 1

/**
* @brief Handles loading input s_data from global memory to shared memory 
* (vec4 version)
*
* Load a chunk of 8*blockDim.x elements from global memory into a 
* shared memory array.  Each thread loads two T4 elements (where
* T4 is, e.g. int4 or float4), computes the scan of those two vec4s in 
* thread local arrays (in registers), and writes the two total sums of the
* vec4s into shared memory, where they will be cooperatively scanned with 
* the other partial sums by all threads in the CTA.
*
* @param[out] s_out The output (shared) memory array
* @param[out] threadScan0 Intermediate per-thread partial sums array 1
* @param[out] threadScan1 Intermediate per-thread partial sums array 2
* @param[in] d_in The input (device) memory array
* @param[in] numElements The number of elements in the array being scanned
* @param[in] iDataOffset the offset of the input array in global memory for this 
* thread block
* @param[out] ai The shared memory address for the thread's first element 
* (returned for reuse)
* @param[out] bi The shared memory address for the thread's second element 
* (returned for reuse)
* @param[out] aiDev The device memory address for this thread's first element 
* (returned for reuse)
* @param[out] biDev The device memory address for this thread's second element 
* (returned for reuse)
*/
template <class T, class traits> 
__device__ void loadSharedChunkFromMem4(T        *s_out,
                                        T        threadScan0[4],
                                        T        threadScan1[4],
                                        const T  *d_in,
                                        int      numElements, 
                                        int      iDataOffset,
                                        int      &ai, 
                                        int      &bi, 
                                        int      &aiDev, 
                                        int      &biDev)
{
    int thid = threadIdx.x;
    aiDev = iDataOffset + thid;
    biDev = aiDev + blockDim.x;

    // convert to 4-vector
    typename typeToVector<T,4>::Result  tempData;
    typename typeToVector<T,4>::Result* inData = (typename typeToVector<T,4>::Result*)d_in;

    ai = thid;
    bi = thid + blockDim.x;

    // read into tempData;
    if (traits::isBackward())
    {
        int i = aiDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements) 
        {
            tempData       = inData[aiDev];
            threadScan0[3] = tempData.w;               
            threadScan0[2] = traits::op(tempData.z, threadScan0[3]);
            threadScan0[1] = traits::op(tempData.y, threadScan0[2]);
            threadScan0[0] = s_out[ai] 
                           = traits::op(tempData.x, threadScan0[1]);
        }
        else
        {
            threadScan0[3] = traits::identity();
            threadScan0[2] = traits::op(((i+2) < numElements) ? d_in[i+2] : traits::identity(), threadScan0[3]);
            threadScan0[1] = traits::op(((i+1) < numElements) ? d_in[i+1] : traits::identity(), threadScan0[2]);
            threadScan0[0] = s_out[ai] 
                           = traits::op((i     < numElements) ? d_in[i]   : traits::identity(), threadScan0[1]);
        }

#ifdef DISALLOW_LOADSTORE_OVERLAP
        __syncthreads();
#endif

        i = biDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements)
        {
            tempData       = inData[biDev];
            threadScan1[3] = tempData.w;
            threadScan1[2] = traits::op(tempData.z, threadScan1[3]);
            threadScan1[1] = traits::op(tempData.y, threadScan1[2]);
            threadScan1[0] = s_out[bi] 
                           = traits::op(tempData.x, threadScan1[1]);
        }
        else
        {
            threadScan1[3] = traits::identity();
            threadScan1[2] = traits::op(((i+2) < numElements) ? d_in[i+2] : traits::identity(), threadScan1[3]);
            threadScan1[1] = traits::op(((i+1) < numElements) ? d_in[i+1] : traits::identity(), threadScan1[2]);
            threadScan1[0] = s_out[bi] 
                           = traits::op((i     < numElements) ? d_in[i]   : traits::identity(), threadScan1[1]);
        }
        __syncthreads();

        // reverse s_data in shared memory
        if (ai < CTA_SIZE)
        {       
            unsigned int leftIdx = ai;
            unsigned int rightIdx = (2 * CTA_SIZE - 1) - ai;
                
            if (leftIdx < rightIdx) 
            {
                T tmp           = s_out[leftIdx];
                s_out[leftIdx]  = s_out[rightIdx];
                s_out[rightIdx] = tmp;
            }
        }
        __syncthreads();
    }
    else
    {
        int i = aiDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements)
        {
            tempData       = inData[aiDev];
            threadScan0[0] = tempData.x;           
            threadScan0[1] = traits::op(tempData.y, threadScan0[0]);
            threadScan0[2] = traits::op(tempData.z, threadScan0[1]);
            threadScan0[3] = s_out[ai] 
                           = traits::op(tempData.w, threadScan0[2]);
        }
        else
        {
            threadScan0[0] = (i < numElements) ? d_in[i] : traits::identity();
            threadScan0[1] = traits::op(((i+1) < numElements) ? d_in[i+1] : traits::identity(), threadScan0[0]);
            threadScan0[2] = traits::op(((i+2) < numElements) ? d_in[i+2] : traits::identity(), threadScan0[1]);
            threadScan0[3] = s_out[ai] 
                           = traits::op(((i+3) < numElements) ? d_in[i+3] : traits::identity(), threadScan0[2]);
        }

        
#ifdef DISALLOW_LOADSTORE_OVERLAP
        __syncthreads();
#endif

        i = biDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements)
        {
            tempData       = inData[biDev];
            threadScan1[0] = tempData.x;           
            threadScan1[1] = traits::op(tempData.y, threadScan1[0]);
            threadScan1[2] = traits::op(tempData.z, threadScan1[1]);
            threadScan1[3] = s_out[bi] 
                           = traits::op(tempData.w, threadScan1[2]);
        }
        else
        {
            threadScan1[0] = (i < numElements) ? d_in[i] : traits::identity();
            threadScan1[1] = traits::op(((i+1) < numElements) ? d_in[i+1] : traits::identity(), threadScan1[0]);
            threadScan1[2] = traits::op(((i+2) < numElements) ? d_in[i+2] : traits::identity(), threadScan1[1]);
            threadScan1[3] = s_out[bi] 
                           = traits::op(((i+3) < numElements) ? d_in[i+3] : traits::identity(), threadScan1[2]);
        }  
        __syncthreads();
    }
}


/**
* @brief Handles storing result s_data from shared memory to global memory 
* (vec4 version)
*
* Store a chunk of SCAN_ELTS_PER_THREAD*blockDim.x elements from shared memory 
* into a device memory array.  Each thread stores reads two elements from shared
* memory, adds them to the intermediate sums computed in 
* loadSharedChunkFromMem4(), and writes two T4 elements (where
* T4 is, e.g. int4 or float4) to global memory.
*
* @param[out] d_out The output (device) memory array
* @param[in] threadScan0 Intermediate per-thread partial sums array 1
* (contents computed in loadSharedChunkFromMem4())
* @param[in] threadScan1 Intermediate per-thread partial sums array 2
* (contents computed in loadSharedChunkFromMem4())
* @param[in] s_in The input (shared) memory array
* @param[in] numElements The number of elements in the array being scanned
* @param[in] oDataOffset the offset of the output array in global memory 
* for this thread block
* @param[in] ai The shared memory address for the thread's first element 
* (computed in loadSharedChunkFromMem4())
* @param[in] bi The shared memory address for the thread's second element 
* (computed in loadSharedChunkFromMem4())
* @param[in] aiDev The device memory address for this thread's first element 
* (computed in loadSharedChunkFromMem4())
* @param[in] biDev The device memory address for this thread's second element 
* (computed in loadSharedChunkFromMem4())
*/
template <class T, class traits>
__device__ void storeSharedChunkToMem4(T   *d_out,
                                       T   threadScan0[4],
                                       T   threadScan1[4],
                                       T   *s_in,
                                       int numElements, 
                                       int oDataOffset,
                                       int ai, 
                                       int bi, 
                                       int aiDev, 
                                       int biDev)
{
    // Convert to 4-vector
    typename typeToVector<T,4>::Result tempData;
    typename typeToVector<T,4>::Result* outData = (typename typeToVector<T,4>::Result*)d_out;

    // write results to global memory
    if (traits::isBackward())
    {   
        if (ai < CTA_SIZE)
        {

            unsigned int leftIdx = ai;
            unsigned int rightIdx = (2 * CTA_SIZE - 1) - ai;
            
            if (leftIdx < rightIdx) 
            {
                T tmp = s_in[leftIdx];
                s_in[leftIdx] = s_in[rightIdx];
                s_in[rightIdx] = tmp;
            }
        }
        __syncthreads();

        T temp = s_in[ai];

        if (traits::isExclusive())
        {
            tempData.w = temp;
            tempData.z = traits::op(temp, threadScan0[3]);
            tempData.y = traits::op(temp, threadScan0[2]);
            tempData.x = traits::op(temp, threadScan0[1]);
        }
        else
        {
            tempData.w = traits::op(temp, threadScan0[3]);
            tempData.z = traits::op(temp, threadScan0[2]);
            tempData.y = traits::op(temp, threadScan0[1]);
            tempData.x = traits::op(temp, threadScan0[0]);
        }

        int i = aiDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements)
        {
            outData[aiDev] = tempData;
        }
        else
        {
            if (i   < numElements) { d_out[i]   = tempData.x;
            if (i+1 < numElements) { d_out[i+1] = tempData.y;
            if (i+2 < numElements) { d_out[i+2] = tempData.z; }}}     
        }

#ifdef DISALLOW_LOADSTORE_OVERLAP
        __syncthreads();
#endif

        temp = s_in[bi];

        if (traits::isExclusive())
        {
            tempData.w = temp;
            tempData.z = traits::op(temp, threadScan1[3]);
            tempData.y = traits::op(temp, threadScan1[2]);
            tempData.x = traits::op(temp, threadScan1[1]);
        }
        else
        {
            tempData.w = traits::op(temp, threadScan1[3]);
            tempData.z = traits::op(temp, threadScan1[2]);
            tempData.y = traits::op(temp, threadScan1[1]);
            tempData.x = traits::op(temp, threadScan1[0]);
        }

        i = biDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements)
        {
            outData[biDev] = tempData;
        }
        else
        {
            if (i   < numElements) { d_out[i]   = tempData.x;
            if (i+1 < numElements) { d_out[i+1] = tempData.y;
            if (i+2 < numElements) { d_out[i+2] = tempData.z; }}}     
        }
    }
    else
    {
        T temp;
        temp = s_in[ai]; 

        if (traits::isExclusive())
        {
            tempData.x = temp;
            tempData.y = traits::op(temp, threadScan0[0]);
            tempData.z = traits::op(temp, threadScan0[1]);
            tempData.w = traits::op(temp, threadScan0[2]);
        }
        else
        {
            tempData.x = traits::op(temp, threadScan0[0]);
            tempData.y = traits::op(temp, threadScan0[1]);
            tempData.z = traits::op(temp, threadScan0[2]);
            tempData.w = traits::op(temp, threadScan0[3]);
        }

        int i = aiDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements)
        {                       
            outData[aiDev] = tempData; 
        }
        else 
        {       
            // we can't use vec4 because the original array isn't a multiple of 
            // 4 elements
            if ( i    < numElements) { d_out[i]   = tempData.x;
            if ((i+1) < numElements) { d_out[i+1] = tempData.y;
            if ((i+2) < numElements) { d_out[i+2] = tempData.z; } } }
        }

#ifdef DISALLOW_LOADSTORE_OVERLAP
        __syncthreads();
#endif

        temp       = s_in[bi]; 

        if (traits::isExclusive())
        {
            tempData.x = temp;
            tempData.y = traits::op(temp, threadScan1[0]);
            tempData.z = traits::op(temp, threadScan1[1]);
            tempData.w = traits::op(temp, threadScan1[2]);
        }
        else
        {
            tempData.x = traits::op(temp, threadScan1[0]);
            tempData.y = traits::op(temp, threadScan1[1]);
            tempData.z = traits::op(temp, threadScan1[2]);
            tempData.w = traits::op(temp, threadScan1[3]);
        }

        i = biDev * 4;
        if (traits::isFullBlock() || i + 3 < numElements)
        {
            outData[biDev] = tempData;
        }
        else 
        {
            // we can't use vec4 because the original array isn't a multiple of 
            // 4 elements
            if ( i    < numElements) { d_out[i]   = tempData.x;
            if ((i+1) < numElements) { d_out[i+1] = tempData.y;
            if ((i+2) < numElements) { d_out[i+2] = tempData.z; } } }
        }
    }
}

/** @brief Scan all warps of a CTA without synchronization
  * 
  * The warp-scan algorithm breaks a block of data into warp-sized chunks, and
  * scans the chunks independently with a warp of threads each.  Because warps
  * execute instructions in SIMD fashion, there is no need to synchronize in 
  * order to share data within a warp (only across warps).  Also, in SIMD the 
  * most efficient algorithm is a step-efficient algorithm.  Therefore, within
  * each warp we use a Hillis-and-Steele-style scan that takes log2(N) steps
  * to scan the warp [Daniel Hillis and Guy Steele 1986], rather than the 
  * work-efficient tree-based algorithm described by Guy Blelloch [1990] that 
  * takes 2 * log(N) steps and is in general more complex to implement.  
  * Previous versions of CUDPP used the Blelloch algorithm.  For current GPUs, 
  * the warp size is 32, so this takes five steps per warp.
  *
  * Each thread is responsible for a single element of the array to be scanned.
  * Each thread inputs a single value to the scan via \a val and returns 
  * its own scanned result element.  The threads of each warp cooperate 
  * via the shared memory array \a s_data to scan WARP_SIZE elements.
  *
  * Template parameter \a maxlevel allows this warpscan to be performed on
  * partial warps.  For example, if only the first 8 elements of each warp need
  * to be scanned, then warpscan only performs log2(8)=3 steps rather than 5.
  *
  * The computation uses 2 * WARP_SIZE elements of shared memory per warp to
  * enable warps to offset beyond their input data and receive the identity 
  * element without using any branch instructions.
  * 
  * \note s_data is declared volatile here to prevent the compiler from 
  * optimizing away writes to shared memory, and ensure correct intrawarp 
  * communication in the absence of __syncthreads.
  *
  * @return The result of the warp scan for the current thread
  * @param[in] val The current threads's input to the scan
  * @param[in,out] s_data A pointer to a temporary shared array of 2*CTA_SIZE 
  * elements used to compute the warp scans
  */
template<class T, class traits,int maxlevel>
__device__ T warpscan(T val, volatile T* s_data)
{
    // The following is the same as 2 * 32 * warpId + threadInWarp = 
    // 64*(threadIdx.x >> 5) + (threadIdx.x & (WARP_SIZE-1))
    int idx = 2 * threadIdx.x - (threadIdx.x & (WARP_SIZE-1));
    s_data[idx] = traits::identity();
    idx += WARP_SIZE;
    T t = s_data[idx] = val;  __EMUSYNC;

        // This code is needed because the warp size of device emulation
        // is only 1 thread, so sync-less cooperation within a warp doesn't 
        // work.
#ifdef __DEVICE_EMULATION__
    t = s_data[idx -  1]; __EMUSYNC; 
    s_data[idx] = traits::op(s_data[idx],t); __EMUSYNC;
    t = s_data[idx -  2]; __EMUSYNC; 
    s_data[idx] = traits::op(s_data[idx],t); __EMUSYNC;
    t = s_data[idx -  4]; __EMUSYNC; 
    s_data[idx] = traits::op(s_data[idx],t); __EMUSYNC;
    t = s_data[idx -  8]; __EMUSYNC; 
    s_data[idx] = traits::op(s_data[idx],t); __EMUSYNC;
    t = s_data[idx - 16]; __EMUSYNC; 
    s_data[idx] = traits::op(s_data[idx],t); __EMUSYNC;
#else
    if (0 <= maxlevel) { s_data[idx] = t = traits::op(t, s_data[idx - 1]); }
    if (1 <= maxlevel) { s_data[idx] = t = traits::op(t, s_data[idx - 2]); }
    if (2 <= maxlevel) { s_data[idx] = t = traits::op(t, s_data[idx - 4]); }
    if (3 <= maxlevel) { s_data[idx] = t = traits::op(t, s_data[idx - 8]); }
    if (4 <= maxlevel) { s_data[idx] = t = traits::op(t, s_data[idx -16]); }
#endif

    return s_data[idx-1];      // convert inclusive -> exclusive
}

/** @brief Perform a full CTA scan using the warp-scan algorithm
  * 
  * As described in the comment for warpscan(), the warp-scan algorithm breaks 
  * a block of data into warp-sized chunks, and scans the chunks independently 
  * with a warp of threads each.  To complete the scan, each warp <i>j</i> then 
  * writes its last element to element <i>j</i> of a temporary shared array.
  * Then a single warp exclusive-scans these "warp sums".  Finally, each thread
  * adds the result of the warp sum scan to the result of the scan from the 
  * first pass.
  *
  * Because we scan 2*CTA_SIZE elements per thread, we have to call warpscan
  * twice.
  *
  * @param x The first input value for the current thread
  * @param y The second input value for the current thread
  * @param s_data Temporary shared memory space of 2*CTA_SIZE elements for 
  * performing the scan
  */
template <class T, class traits>
__device__ void scanWarps(T x, T y, 
                          T *s_data)
{       
    T val  = warpscan<T, traits, 4>(x, s_data);
    __syncthreads(); 
    T val2 = warpscan<T, traits, 4>(y, s_data);
    
    int idx = threadIdx.x;

    if ((idx & 31)==31)
    {
        s_data[idx >> 5]                = traits::op(val, x);
        s_data[(idx + blockDim.x) >> 5] = traits::op(val2, y);
    }
    __syncthreads();

#ifndef __DEVICE_EMULATION__
    if (idx < 32)
#endif
    {
        s_data[idx] = warpscan<T,traits,(LOG_CTA_SIZE-LOG_WARP_SIZE+1)>(s_data[idx], s_data);
    }
    __syncthreads();

    val  = traits::op(val, s_data[idx >> 5]);

    val2 = traits::op(val2, s_data[(idx + blockDim.x) >> 5]);

    __syncthreads();

    s_data[idx] = val;
    s_data[idx+blockDim.x] = val2;
}

/**
* @brief CTA-level scan routine; scans s_data in shared memory in each thread block
*
* This function is the main CTA-level scan function.  It may be called by other 
* CUDA __global__ or __device__ functions. This function scans 2 * CTA_SIZE elements.
* Each thread is responsible for one element in each half of the input array.
* \note This code is intended to be run on a CTA of 128 threads.  Other sizes are
* untested.
* 
* @param[in] s_data The array to be scanned in shared memory
* @param[out] d_blockSums Array of per-block sums
* @param[in] blockSumIndex Location in \a d_blockSums to which to write this block's sum
*/
template <class T, class traits>
__device__ void scanCTA(T            *s_data, 
                        T            *d_blockSums, 
                        unsigned int blockSumIndex)
{
    T val  = s_data[threadIdx.x];
    T val2 = s_data[threadIdx.x + blockDim.x];
    __syncthreads();     

    scanWarps<T,traits>(val, val2, s_data);
    __syncthreads();  

    if (traits::writeSums() && threadIdx.x == blockDim.x - 1)
    {
        d_blockSums[blockSumIndex] = traits::op(val2, s_data[threadIdx.x + blockDim.x]);
    }
    
    
#ifdef __DEVICE_EMULATION__
    // must sync in emulation mode when doing backward scans, because otherwise the 
    // shared memory array will get reversed before the block sums are read!
    if (traits::isBackward())
        __syncthreads();
#endif
}


/** @} */ // end scan functions
/** @} */ // end cudpp_cta
