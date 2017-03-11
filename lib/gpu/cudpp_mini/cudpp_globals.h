// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 5289 $
// $Date: 2010-11-23 13:04:43 -0700 (Tue, 23 Nov 2010) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt in
// the root directory of this source distribution.
// ------------------------------------------------------------- 
 
/**
 * @file
 * cudpp_globals.h
 *
 * @brief Global declarations defining machine characteristics of GPU target
 * These are currently set for best performance on G8X GPUs.  The optimal 
 * parameters may change on future GPUs. In the future, we hope to make
 * CUDPP a self-tuning library.
 */

#ifndef __CUDPP_GLOBALS_H__
#define __CUDPP_GLOBALS_H__

const int NUM_BANKS = 16;                        /**< Number of shared memory banks */
const int LOG_NUM_BANKS = 4;                     /**< log_2(NUM_BANKS) */
const int CTA_SIZE = 128;                        /**< Number of threads in a CTA */
const int WARP_SIZE = 32;                        /**< Number of threads in a warp */
const int LOG_CTA_SIZE = 7;                      /**< log_2(CTA_SIZE) */
const int LOG_WARP_SIZE = 5;                     /**< log_2(WARP_SIZE) */
const int LOG_SIZEOF_FLOAT = 2;                  /**< log_2(sizeof(float)) */
const int SCAN_ELTS_PER_THREAD = 8;              /**< Number of elements per scan thread */
const int SEGSCAN_ELTS_PER_THREAD = 8;     /**< Number of elements per segmented scan thread */

const int maxSharedMemoryPerBlock = 16384; /**< Number of bytes of shared 
                                              memory in each block */
const int maxThreadsPerBlock = CTA_SIZE;   /**< Maximum number of
                                             * threads in a CTA */

/**
* @brief Macro to insert necessary __syncthreads() in device emulation mode
*/
#ifdef __DEVICE_EMULATION__
#define __EMUSYNC __syncthreads()
#else
#define __EMUSYNC
#endif


#define AVOID_BANK_CONFLICTS /**< Set if by default, we want our
                              * shared memory allocation to perform
                              * additional computation to avoid bank
                              * conflicts */

#ifdef AVOID_BANK_CONFLICTS
#define CONFLICT_FREE_OFFSET(index) ((index) >> LOG_NUM_BANKS)
#else
#define CONFLICT_FREE_OFFSET(index) (0)
#endif

#endif // __CUDPP_GLOBALS_H__

// Leave this at the end of the file
// Local Variables:
// mode:c++
// c-file-style: "NVIDIA"
// End:
