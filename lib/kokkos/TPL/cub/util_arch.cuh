/******************************************************************************
 * Copyright (c) 2011, Duane Merrill.  All rights reserved.
 * Copyright (c) 2011-2013, NVIDIA CORPORATION.  All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

/**
 * \file
 * Static architectural properties by SM version.
 */


/******************************************************************************
 * Static architectural properties by SM version.
 *
 * "Device" reflects the PTX architecture targeted by the active compiler
 * pass.  It provides useful compile-time statics within device code.  E.g.,:
 *
 *     __shared__ int[Device::WARP_THREADS];
 *
 *     int padded_offset = threadIdx.x + (threadIdx.x >> Device::LOG_SMEM_BANKS);
 *
 ******************************************************************************/

#pragma once

#include "util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/**
 * \addtogroup UtilModule
 * @{
 */


/// CUB_PTX_ARCH reflects the PTX version targeted by the active compiler pass (or zero during the host pass).
#ifndef __CUDA_ARCH__
    #define CUB_PTX_ARCH 0
#else
    #define CUB_PTX_ARCH __CUDA_ARCH__
#endif


/// Whether or not the source targeted by the active compiler pass is allowed to  invoke device kernels or methods from the CUDA runtime API.
#if !defined(__CUDA_ARCH__) || defined(CUB_CDP)
#define CUB_RUNTIME_ENABLED
#endif


/// Execution space for destructors
#if ((CUB_PTX_ARCH > 0) && (CUB_PTX_ARCH < 200))
    #define CUB_DESTRUCTOR __host__
#else
    #define CUB_DESTRUCTOR __host__ __device__
#endif


/**
 * \brief Structure for statically reporting CUDA device properties, parameterized by SM architecture.
 *
 * The default specialization is for SM10.
 */
template <int SM_ARCH>
struct ArchProps
{
    enum
    {
        LOG_WARP_THREADS    =
                                        5,                        /// Log of the number of threads per warp
        WARP_THREADS        =
                                        1 << LOG_WARP_THREADS,    /// Number of threads per warp
        LOG_SMEM_BANKS      =
                                        4,                        /// Log of the number of smem banks
        SMEM_BANKS          =
                                        1 << LOG_SMEM_BANKS,      /// The number of smem banks
        SMEM_BANK_BYTES     =
                                        4,                        /// Size of smem bank words
        SMEM_BYTES          =
                                        16 * 1024,                /// Maximum SM shared memory
        SMEM_ALLOC_UNIT     =
                                        512,                      /// Smem allocation size in bytes
        REGS_BY_BLOCK       =
                                        true,                     /// Whether or not the architecture allocates registers by block (or by warp)
        REG_ALLOC_UNIT      =
                                        256,                      /// Number of registers allocated at a time per block (or by warp)
        WARP_ALLOC_UNIT     =
                                        2,                        /// Granularity of warps for which registers are allocated
        MAX_SM_THREADS      =
                                        768,                      /// Maximum number of threads per SM
        MAX_SM_THREADBLOCKS =
                                        8,                        /// Maximum number of thread blocks per SM
        MAX_BLOCK_THREADS   =
                                        512,                      /// Maximum number of thread per thread block
        MAX_SM_REGISTERS    =
                                        8 * 1024,                 /// Maximum number of registers per SM
    };
};




#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

/**
 * Architecture properties for SM30
 */
template <>
struct ArchProps<300>
{
    enum
    {
        LOG_WARP_THREADS    = 5,                        // 32 threads per warp
        WARP_THREADS        = 1 << LOG_WARP_THREADS,
        LOG_SMEM_BANKS      = 5,                        // 32 banks
        SMEM_BANKS          = 1 << LOG_SMEM_BANKS,
        SMEM_BANK_BYTES     = 4,                        // 4 byte bank words
        SMEM_BYTES          = 48 * 1024,                // 48KB shared memory
        SMEM_ALLOC_UNIT     = 256,                      // 256B smem allocation segment size
        REGS_BY_BLOCK       = false,                    // Allocates registers by warp
        REG_ALLOC_UNIT      = 256,                      // 256 registers allocated at a time per warp
        WARP_ALLOC_UNIT     = 4,                        // Registers are allocated at a granularity of every 4 warps per threadblock
        MAX_SM_THREADS      = 2048,                     // 2K max threads per SM
        MAX_SM_THREADBLOCKS = 16,                       // 16 max threadblocks per SM
        MAX_BLOCK_THREADS   = 1024,                     // 1024 max threads per threadblock
        MAX_SM_REGISTERS    = 64 * 1024,                // 64K max registers per SM
    };

    // Callback utility
    template <typename T>
    static __host__ __device__ __forceinline__ void Callback(T &target, int sm_version)
    {
        target.template Callback<ArchProps>();
    }
};


/**
 * Architecture properties for SM20
 */
template <>
struct ArchProps<200>
{
    enum
    {
        LOG_WARP_THREADS    = 5,                        // 32 threads per warp
        WARP_THREADS        = 1 << LOG_WARP_THREADS,
        LOG_SMEM_BANKS      = 5,                        // 32 banks
        SMEM_BANKS          = 1 << LOG_SMEM_BANKS,
        SMEM_BANK_BYTES     = 4,                        // 4 byte bank words
        SMEM_BYTES          = 48 * 1024,                // 48KB shared memory
        SMEM_ALLOC_UNIT     = 128,                      // 128B smem allocation segment size
        REGS_BY_BLOCK       = false,                    // Allocates registers by warp
        REG_ALLOC_UNIT      = 64,                       // 64 registers allocated at a time per warp
        WARP_ALLOC_UNIT     = 2,                        // Registers are allocated at a granularity of every 2 warps per threadblock
        MAX_SM_THREADS      = 1536,                     // 1536 max threads per SM
        MAX_SM_THREADBLOCKS = 8,                        // 8 max threadblocks per SM
        MAX_BLOCK_THREADS   = 1024,                     // 1024 max threads per threadblock
        MAX_SM_REGISTERS    = 32 * 1024,                // 32K max registers per SM
    };

    // Callback utility
    template <typename T>
    static __host__ __device__ __forceinline__ void Callback(T &target, int sm_version)
    {
        if (sm_version > 200) {
            ArchProps<300>::Callback(target, sm_version);
        } else {
            target.template Callback<ArchProps>();
        }
    }
};


/**
 * Architecture properties for SM12
 */
template <>
struct ArchProps<120>
{
    enum
    {
        LOG_WARP_THREADS    = 5,                        // 32 threads per warp
        WARP_THREADS        = 1 << LOG_WARP_THREADS,
        LOG_SMEM_BANKS      = 4,                        // 16 banks
        SMEM_BANKS          = 1 << LOG_SMEM_BANKS,
        SMEM_BANK_BYTES     = 4,                        // 4 byte bank words
        SMEM_BYTES          = 16 * 1024,                // 16KB shared memory
        SMEM_ALLOC_UNIT     = 512,                      // 512B smem allocation segment size
        REGS_BY_BLOCK       = true,                     // Allocates registers by threadblock
        REG_ALLOC_UNIT      = 512,                      // 512 registers allocated at time per threadblock
        WARP_ALLOC_UNIT     = 2,                        // Registers are allocated at a granularity of every 2 warps per threadblock
        MAX_SM_THREADS      = 1024,                     // 1024 max threads per SM
        MAX_SM_THREADBLOCKS = 8,                        // 8 max threadblocks per SM
        MAX_BLOCK_THREADS   = 512,                      // 512 max threads per threadblock
        MAX_SM_REGISTERS    = 16 * 1024,                // 16K max registers per SM
    };

    // Callback utility
    template <typename T>
    static __host__ __device__ __forceinline__ void Callback(T &target, int sm_version)
    {
        if (sm_version > 120) {
            ArchProps<200>::Callback(target, sm_version);
        } else {
            target.template Callback<ArchProps>();
        }
    }
};


/**
 * Architecture properties for SM10.  Derives from the default ArchProps specialization.
 */
template <>
struct ArchProps<100> : ArchProps<0>
{
    // Callback utility
    template <typename T>
    static __host__ __device__ __forceinline__ void Callback(T &target, int sm_version)
    {
        if (sm_version > 100) {
            ArchProps<120>::Callback(target, sm_version);
        } else {
            target.template Callback<ArchProps>();
        }
    }
};


/**
 * Architecture properties for SM35
 */
template <>
struct ArchProps<350> : ArchProps<300> {};        // Derives from SM30

/**
 * Architecture properties for SM21
 */
template <>
struct ArchProps<210> : ArchProps<200> {};        // Derives from SM20

/**
 * Architecture properties for SM13
 */
template <>
struct ArchProps<130> : ArchProps<120> {};        // Derives from SM12

/**
 * Architecture properties for SM11
 */
template <>
struct ArchProps<110> : ArchProps<100> {};        // Derives from SM10


#endif // DOXYGEN_SHOULD_SKIP_THIS


/**
 * \brief The architectural properties for the PTX version targeted by the active compiler pass.
 */
struct PtxArchProps : ArchProps<CUB_PTX_ARCH> {};


/** @} */       // end group UtilModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
