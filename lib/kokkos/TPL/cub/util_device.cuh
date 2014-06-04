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
 * Properties of a given CUDA device and the corresponding PTX bundle
 */

#pragma once

#include "util_arch.cuh"
#include "util_debug.cuh"
#include "util_namespace.cuh"
#include "util_macro.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/**
 * \addtogroup UtilModule
 * @{
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document


/**
 * Empty kernel for querying PTX manifest metadata (e.g., version) for the current device
 */
template <typename T>
__global__ void EmptyKernel(void) { }


/**
 * Alias temporaries to externally-allocated device storage (or simply return the amount of storage needed).
 */
template <int ALLOCATIONS>
__host__ __device__ __forceinline__
cudaError_t AliasTemporaries(
    void    *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
    size_t  &temp_storage_bytes,                ///< [in,out] Size in bytes of \t d_temp_storage allocation
    void*   (&allocations)[ALLOCATIONS],        ///< [in,out] Pointers to device allocations needed
    size_t  (&allocation_sizes)[ALLOCATIONS])   ///< [in] Sizes in bytes of device allocations needed
{
    const int ALIGN_BYTES   = 256;
    const int ALIGN_MASK    = ~(ALIGN_BYTES - 1);

    // Compute exclusive prefix sum over allocation requests
    size_t bytes_needed = 0;
    for (int i = 0; i < ALLOCATIONS; ++i)
    {
        size_t allocation_bytes = (allocation_sizes[i] + ALIGN_BYTES - 1) & ALIGN_MASK;
        allocation_sizes[i] = bytes_needed;
        bytes_needed += allocation_bytes;
    }

    // Check if the caller is simply requesting the size of the storage allocation
    if (!d_temp_storage)
    {
        temp_storage_bytes = bytes_needed;
        return cudaSuccess;
    }

    // Check if enough storage provided
    if (temp_storage_bytes < bytes_needed)
    {
        return CubDebug(cudaErrorMemoryAllocation);
    }

    // Alias
    for (int i = 0; i < ALLOCATIONS; ++i)
    {
        allocations[i] = static_cast<char*>(d_temp_storage) + allocation_sizes[i];
    }

    return cudaSuccess;
}



#endif  // DOXYGEN_SHOULD_SKIP_THIS



/**
 * \brief Retrieves the PTX version (major * 100 + minor * 10)
 */
__host__ __device__ __forceinline__ cudaError_t PtxVersion(int &ptx_version)
{
#ifndef CUB_RUNTIME_ENABLED

    // CUDA API calls not supported from this device
    return cudaErrorInvalidConfiguration;

#else

    cudaError_t error = cudaSuccess;
    do
    {
        cudaFuncAttributes empty_kernel_attrs;
        if (CubDebug(error = cudaFuncGetAttributes(&empty_kernel_attrs, EmptyKernel<void>))) break;
        ptx_version = empty_kernel_attrs.ptxVersion * 10;
    }
    while (0);

    return error;

#endif
}


/**
 * Synchronize the stream if specified
 */
__host__ __device__ __forceinline__
static cudaError_t SyncStream(cudaStream_t stream)
{
#ifndef __CUDA_ARCH__
    return cudaStreamSynchronize(stream);
#else
    // Device can't yet sync on a specific stream
    return cudaDeviceSynchronize();
#endif
}



/**
 * \brief Properties of a given CUDA device and the corresponding PTX bundle
 */
class Device
{
private:

    /// Type definition of the EmptyKernel kernel entry point
    typedef void (*EmptyKernelPtr)();

    /// Force EmptyKernel<void> to be generated if this class is used
    __host__ __device__ __forceinline__
    EmptyKernelPtr Empty()
    {
        return EmptyKernel<void>;
    }

public:

    // Version information
    int     sm_version;             ///< SM version of target device (SM version X.YZ in XYZ integer form)
    int     ptx_version;            ///< Bundled PTX version for target device (PTX version X.YZ in XYZ integer form)

    // Target device properties
    int     sm_count;               ///< Number of SMs
    int     warp_threads;           ///< Number of threads per warp
    int     smem_bank_bytes;        ///< Number of bytes per SM bank
    int     smem_banks;             ///< Number of smem banks
    int     smem_bytes;             ///< Smem bytes per SM
    int     smem_alloc_unit;        ///< Smem segment size
    bool    regs_by_block;          ///< Whether registers are allocated by threadblock (or by warp)
    int     reg_alloc_unit;         ///< Granularity of register allocation within the SM
    int     warp_alloc_unit;        ///< Granularity of warp allocation within the SM
    int     max_sm_threads;         ///< Maximum number of threads per SM
    int     max_sm_blocks;          ///< Maximum number of threadblocks per SM
    int     max_block_threads;      ///< Maximum number of threads per threadblock
    int     max_sm_registers;       ///< Maximum number of registers per SM
    int     max_sm_warps;           ///< Maximum number of warps per SM

    /**
     * Callback for initializing device properties
     */
    template <typename ArchProps>
    __host__ __device__ __forceinline__ void Callback()
    {
        warp_threads        = ArchProps::WARP_THREADS;
        smem_bank_bytes     = ArchProps::SMEM_BANK_BYTES;
        smem_banks          = ArchProps::SMEM_BANKS;
        smem_bytes          = ArchProps::SMEM_BYTES;
        smem_alloc_unit     = ArchProps::SMEM_ALLOC_UNIT;
        regs_by_block       = ArchProps::REGS_BY_BLOCK;
        reg_alloc_unit      = ArchProps::REG_ALLOC_UNIT;
        warp_alloc_unit     = ArchProps::WARP_ALLOC_UNIT;
        max_sm_threads      = ArchProps::MAX_SM_THREADS;
        max_sm_blocks       = ArchProps::MAX_SM_THREADBLOCKS;
        max_block_threads   = ArchProps::MAX_BLOCK_THREADS;
        max_sm_registers    = ArchProps::MAX_SM_REGISTERS;
        max_sm_warps        = max_sm_threads / warp_threads;
    }


public:

    /**
     * Initializer.  Properties are retrieved for the specified GPU ordinal.
     */
    __host__ __device__ __forceinline__
    cudaError_t Init(int device_ordinal)
    {
    #ifndef CUB_RUNTIME_ENABLED

        // CUDA API calls not supported from this device
        return CubDebug(cudaErrorInvalidConfiguration);

    #else

        cudaError_t error = cudaSuccess;
        do
        {
            // Fill in SM version
            int major, minor;
            if (CubDebug(error = cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, device_ordinal))) break;
            if (CubDebug(error = cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, device_ordinal))) break;
            sm_version = major * 100 + minor * 10;

            // Fill in static SM properties
            // Initialize our device properties via callback from static device properties
            ArchProps<100>::Callback(*this, sm_version);

            // Fill in SM count
            if (CubDebug(error = cudaDeviceGetAttribute (&sm_count, cudaDevAttrMultiProcessorCount, device_ordinal))) break;

            // Fill in PTX version
        #if CUB_PTX_ARCH > 0
            ptx_version = CUB_PTX_ARCH;
        #else
            if (CubDebug(error = PtxVersion(ptx_version))) break;
        #endif

        }
        while (0);

        return error;

    #endif
    }


    /**
     * Initializer.  Properties are retrieved for the current GPU ordinal.
     */
    __host__ __device__ __forceinline__
    cudaError_t Init()
    {
    #ifndef CUB_RUNTIME_ENABLED

        // CUDA API calls not supported from this device
        return CubDebug(cudaErrorInvalidConfiguration);

    #else

        cudaError_t error = cudaSuccess;
        do
        {
            int device_ordinal;
            if ((error = CubDebug(cudaGetDevice(&device_ordinal)))) break;
            if ((error = Init(device_ordinal))) break;
        }
        while (0);
        return error;

    #endif
    }


    /**
     * Computes maximum SM occupancy in thread blocks for the given kernel
     */
    template <typename KernelPtr>
    __host__ __device__ __forceinline__
    cudaError_t MaxSmOccupancy(
        int                 &max_sm_occupancy,          ///< [out] maximum number of thread blocks that can reside on a single SM
        KernelPtr           kernel_ptr,                 ///< [in] Kernel pointer for which to compute SM occupancy
        int                 block_threads)              ///< [in] Number of threads per thread block
    {
    #ifndef CUB_RUNTIME_ENABLED

        // CUDA API calls not supported from this device
        return CubDebug(cudaErrorInvalidConfiguration);

    #else

        cudaError_t error = cudaSuccess;
        do
        {
            // Get kernel attributes
            cudaFuncAttributes kernel_attrs;
            if (CubDebug(error = cudaFuncGetAttributes(&kernel_attrs, kernel_ptr))) break;

            // Number of warps per threadblock
            int block_warps = (block_threads +  warp_threads - 1) / warp_threads;

            // Max warp occupancy
            int max_warp_occupancy = (block_warps > 0) ?
                max_sm_warps / block_warps :
                max_sm_blocks;

            // Maximum register occupancy
            int max_reg_occupancy;
            if ((block_threads == 0) || (kernel_attrs.numRegs == 0))
            {
                // Prevent divide-by-zero
                max_reg_occupancy = max_sm_blocks;
            }
            else if (regs_by_block)
            {
                // Allocates registers by threadblock
                int block_regs = CUB_ROUND_UP_NEAREST(kernel_attrs.numRegs * warp_threads * block_warps, reg_alloc_unit);
                max_reg_occupancy = max_sm_registers / block_regs;
            }
            else
            {
                // Allocates registers by warp
                int sm_sides                = warp_alloc_unit;
                int sm_registers_per_side   = max_sm_registers / sm_sides;
                int regs_per_warp           = CUB_ROUND_UP_NEAREST(kernel_attrs.numRegs * warp_threads, reg_alloc_unit);
                int warps_per_side          = sm_registers_per_side / regs_per_warp;
                int warps                   = warps_per_side * sm_sides;
                max_reg_occupancy           = warps / block_warps;
            }

            // Shared memory per threadblock
            int block_allocated_smem = CUB_ROUND_UP_NEAREST(
                kernel_attrs.sharedSizeBytes,
                smem_alloc_unit);

            // Max shared memory occupancy
            int max_smem_occupancy = (block_allocated_smem > 0) ?
                (smem_bytes / block_allocated_smem) :
                max_sm_blocks;

            // Max occupancy
            max_sm_occupancy = CUB_MIN(
                CUB_MIN(max_sm_blocks, max_warp_occupancy),
                CUB_MIN(max_smem_occupancy, max_reg_occupancy));

//            printf("max_smem_occupancy(%d), max_warp_occupancy(%d), max_reg_occupancy(%d)", max_smem_occupancy, max_warp_occupancy, max_reg_occupancy);

        } while (0);

        return error;

    #endif
    }

};


/** @} */       // end group UtilModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
