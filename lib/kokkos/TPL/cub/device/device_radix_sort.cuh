
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
 * cub::DeviceRadixSort provides operations for computing a device-wide, parallel reduction across data items residing within global memory.
 */

#pragma once

#include <stdio.h>
#include <iterator>

#include "block/block_radix_sort_upsweep_tiles.cuh"
#include "block/block_radix_sort_downsweep_tiles.cuh"
#include "block/block_scan_tiles.cuh"
#include "../grid/grid_even_share.cuh"
#include "../util_debug.cuh"
#include "../util_device.cuh"
#include "../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document




/******************************************************************************
 * Kernel entry points
 *****************************************************************************/

/**
 * Upsweep pass kernel entry point (multi-block).  Computes privatized digit histograms, one per block.
 */
template <
    typename                BlockRadixSortUpsweepTilesPolicy, ///< Tuning policy for cub::BlockRadixSortUpsweepTiles abstraction
    typename                Key,                            ///< Key type
    typename                SizeT>                          ///< Integer type used for global array indexing
__launch_bounds__ (int(BlockRadixSortUpsweepTilesPolicy::BLOCK_THREADS), 1)
__global__ void RadixSortUpsweepKernel(
    Key                     *d_keys,                        ///< [in] Input keys buffer
    SizeT                   *d_spine,                       ///< [out] Privatized (per block) digit histograms (striped, i.e., 0s counts from each block, then 1s counts from each block, etc.)
    SizeT                   num_items,                      ///< [in] Total number of input data items
    int                     current_bit,                    ///< [in] Bit position of current radix digit
    bool                    use_primary_bit_granularity,    ///< [in] Whether nor not to use the primary policy (or the embedded alternate policy for smaller bit granularity)
    bool                    first_pass,                     ///< [in] Whether this is the first digit pass
    GridEvenShare<SizeT>    even_share)                     ///< [in] Descriptor for how to map an even-share of tiles across thread blocks
{

    // Alternate policy for when fewer bits remain
    typedef typename BlockRadixSortUpsweepTilesPolicy::AltPolicy AltPolicy;

    // Parameterize two versions of BlockRadixSortUpsweepTiles type for the current configuration
    typedef BlockRadixSortUpsweepTiles<BlockRadixSortUpsweepTilesPolicy, Key, SizeT>    BlockRadixSortUpsweepTilesT;          // Primary
    typedef BlockRadixSortUpsweepTiles<AltPolicy, Key, SizeT>                           AltBlockRadixSortUpsweepTilesT;       // Alternate (smaller bit granularity)

    // Shared memory storage
    __shared__ union
    {
        typename BlockRadixSortUpsweepTilesT::TempStorage     pass_storage;
        typename AltBlockRadixSortUpsweepTilesT::TempStorage  alt_pass_storage;
    } temp_storage;

    // Initialize even-share descriptor for this thread block
    even_share.BlockInit();

    // Process input tiles (each of the first RADIX_DIGITS threads will compute a count for that digit)
    if (use_primary_bit_granularity)
    {
        // Primary granularity
        SizeT bin_count;
        BlockRadixSortUpsweepTilesT(temp_storage.pass_storage, d_keys, current_bit).ProcessTiles(
            even_share.block_offset,
            even_share.block_oob,
            bin_count);

        // Write out digit counts (striped)
        if (threadIdx.x < BlockRadixSortUpsweepTilesT::RADIX_DIGITS)
        {
            d_spine[(gridDim.x * threadIdx.x) + blockIdx.x] = bin_count;
        }
    }
    else
    {
        // Alternate granularity
        // Process input tiles (each of the first RADIX_DIGITS threads will compute a count for that digit)
        SizeT bin_count;
        AltBlockRadixSortUpsweepTilesT(temp_storage.alt_pass_storage, d_keys, current_bit).ProcessTiles(
            even_share.block_offset,
            even_share.block_oob,
            bin_count);

        // Write out digit counts (striped)
        if (threadIdx.x < AltBlockRadixSortUpsweepTilesT::RADIX_DIGITS)
        {
            d_spine[(gridDim.x * threadIdx.x) + blockIdx.x] = bin_count;
        }
    }
}


/**
 * Spine scan kernel entry point (single-block).  Computes an exclusive prefix sum over the privatized digit histograms
 */
template <
    typename    BlockScanTilesPolicy,   ///< Tuning policy for cub::BlockScanTiles abstraction
    typename    SizeT>                  ///< Integer type used for global array indexing
__launch_bounds__ (int(BlockScanTilesPolicy::BLOCK_THREADS), 1)
__global__ void RadixSortScanKernel(
    SizeT       *d_spine,               ///< [in,out] Privatized (per block) digit histograms (striped, i.e., 0s counts from each block, then 1s counts from each block, etc.)
    int         num_counts)             ///< [in] Total number of bin-counts
{
    // Parameterize the BlockScanTiles type for the current configuration
    typedef BlockScanTiles<BlockScanTilesPolicy, SizeT*, SizeT*, cub::Sum, SizeT, SizeT> BlockScanTilesT;

    // Shared memory storage
    __shared__ typename BlockScanTilesT::TempStorage temp_storage;

    // Block scan instance
    BlockScanTilesT block_scan(temp_storage, d_spine, d_spine, cub::Sum(), SizeT(0)) ;

    // Process full input tiles
    int block_offset = 0;
    RunningBlockPrefixOp<SizeT> prefix_op;
    prefix_op.running_total = 0;
    while (block_offset < num_counts)
    {
        block_scan.ConsumeTile<true, false>(block_offset, prefix_op);
        block_offset += BlockScanTilesT::TILE_ITEMS;
    }
}


/**
 * Downsweep pass kernel entry point (multi-block).  Scatters keys (and values) into corresponding bins for the current digit place.
 */
template <
    typename                BlockRadixSortDownsweepTilesPolicy,   ///< Tuning policy for cub::BlockRadixSortUpsweepTiles abstraction
    typename                Key,                                ///< Key type
    typename                Value,                              ///< Value type
    typename                SizeT>                              ///< Integer type used for global array indexing
__launch_bounds__ (int(BlockRadixSortDownsweepTilesPolicy::BLOCK_THREADS))
__global__ void RadixSortDownsweepKernel(
    Key                     *d_keys_in,                     ///< [in] Input keys ping buffer
    Key                     *d_keys_out,                    ///< [in] Output keys pong buffer
    Value                   *d_values_in,                   ///< [in] Input values ping buffer
    Value                   *d_values_out,                  ///< [in] Output values pong buffer
    SizeT                   *d_spine,                       ///< [in] Scan of privatized (per block) digit histograms (striped, i.e., 0s counts from each block, then 1s counts from each block, etc.)
    SizeT                   num_items,                      ///< [in] Total number of input data items
    int                     current_bit,                    ///< [in] Bit position of current radix digit
    bool                    use_primary_bit_granularity,    ///< [in] Whether nor not to use the primary policy (or the embedded alternate policy for smaller bit granularity)
    bool                    first_pass,                     ///< [in] Whether this is the first digit pass
    bool                    last_pass,                      ///< [in] Whether this is the last digit pass
    GridEvenShare<SizeT>    even_share)                     ///< [in] Descriptor for how to map an even-share of tiles across thread blocks
{

    // Alternate policy for when fewer bits remain
    typedef typename BlockRadixSortDownsweepTilesPolicy::AltPolicy AltPolicy;

    // Parameterize two versions of BlockRadixSortDownsweepTiles type for the current configuration
    typedef BlockRadixSortDownsweepTiles<BlockRadixSortDownsweepTilesPolicy, Key, Value, SizeT>     BlockRadixSortDownsweepTilesT;
    typedef BlockRadixSortDownsweepTiles<AltPolicy, Key, Value, SizeT>                            AltBlockRadixSortDownsweepTilesT;

    // Shared memory storage
    __shared__ union
    {
        typename BlockRadixSortDownsweepTilesT::TempStorage       pass_storage;
        typename AltBlockRadixSortDownsweepTilesT::TempStorage    alt_pass_storage;

    } temp_storage;

    // Initialize even-share descriptor for this thread block
    even_share.BlockInit();

    if (use_primary_bit_granularity)
    {
        // Process input tiles
        BlockRadixSortDownsweepTilesT(temp_storage.pass_storage, num_items, d_spine, d_keys_in, d_keys_out, d_values_in, d_values_out, current_bit).ProcessTiles(
            even_share.block_offset,
            even_share.block_oob);
    }
    else
    {
        // Process input tiles
        AltBlockRadixSortDownsweepTilesT(temp_storage.alt_pass_storage, num_items, d_spine, d_keys_in, d_keys_out, d_values_in, d_values_out, current_bit).ProcessTiles(
            even_share.block_offset,
            even_share.block_oob);
    }
}


#endif // DOXYGEN_SHOULD_SKIP_THIS





/******************************************************************************
 * DeviceRadixSort
 *****************************************************************************/

/**
 * \brief DeviceRadixSort provides operations for computing a device-wide, parallel radix sort across data items residing within global memory. ![](sorting_logo.png)
 * \ingroup DeviceModule
 *
 * \par Overview
 * The [<em>radix sorting method</em>](http://en.wikipedia.org/wiki/Radix_sort) arranges
 * items into ascending order.  It relies upon a positional representation for
 * keys, i.e., each key is comprised of an ordered sequence of symbols (e.g., digits,
 * characters, etc.) specified from least-significant to most-significant.  For a
 * given input sequence of keys and a set of rules specifying a total ordering
 * of the symbolic alphabet, the radix sorting method produces a lexicographic
 * ordering of those keys.
 *
 * \par
 * DeviceRadixSort can sort all of the built-in C++ numeric primitive types, e.g.:
 * <tt>unsigned char</tt>, \p int, \p double, etc.  Although the direct radix sorting
 * method can only be applied to unsigned integral types, BlockRadixSort
 * is able to sort signed and floating-point types via simple bit-wise transformations
 * that ensure lexicographic key ordering.
 *
 * \par Usage Considerations
 * \cdp_class{DeviceRadixSort}
 *
 * \par Performance
 *
 * \image html lsd_sort_perf.png
 *
 */
struct DeviceRadixSort
{
    #ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document


    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    /// Generic structure for encapsulating dispatch properties codified in block policy.
    struct KernelDispachParams
    {
        int                     block_threads;
        int                     items_per_thread;
        cudaSharedMemConfig     smem_config;
        int                     radix_bits;
        int                     alt_radix_bits;
        int                     subscription_factor;
        int                     tile_size;

        template <typename SortBlockPolicy>
        __host__ __device__ __forceinline__
        void InitUpsweepPolicy(int subscription_factor = 1)
        {
            block_threads               = SortBlockPolicy::BLOCK_THREADS;
            items_per_thread            = SortBlockPolicy::ITEMS_PER_THREAD;
            radix_bits                  = SortBlockPolicy::RADIX_BITS;
            alt_radix_bits              = SortBlockPolicy::AltPolicy::RADIX_BITS;
            smem_config                 = cudaSharedMemBankSizeFourByte;
            this->subscription_factor   = subscription_factor;
            tile_size                   = block_threads * items_per_thread;
        }

        template <typename ScanBlockPolicy>
        __host__ __device__ __forceinline__
        void InitScanPolicy()
        {
            block_threads               = ScanBlockPolicy::BLOCK_THREADS;
            items_per_thread            = ScanBlockPolicy::ITEMS_PER_THREAD;
            radix_bits                  = 0;
            alt_radix_bits              = 0;
            smem_config                 = cudaSharedMemBankSizeFourByte;
            subscription_factor         = 0;
            tile_size                   = block_threads * items_per_thread;
        }

        template <typename SortBlockPolicy>
        __host__ __device__ __forceinline__
        void InitDownsweepPolicy(int subscription_factor = 1)
        {
            block_threads               = SortBlockPolicy::BLOCK_THREADS;
            items_per_thread            = SortBlockPolicy::ITEMS_PER_THREAD;
            radix_bits                  = SortBlockPolicy::RADIX_BITS;
            alt_radix_bits              = SortBlockPolicy::AltPolicy::RADIX_BITS;
            smem_config                 = SortBlockPolicy::SMEM_CONFIG;
            this->subscription_factor   = subscription_factor;
            tile_size                   = block_threads * items_per_thread;
        }
    };



    /******************************************************************************
     * Tuning policies
     ******************************************************************************/

    /// Specializations of tuned policy types for different PTX architectures
    template <typename Key, typename Value, typename SizeT, int ARCH>
    struct TunedPolicies;

    /// SM35 tune
    template <typename Key, typename Value, typename SizeT>
    struct TunedPolicies<Key, Value, SizeT, 350>
    {
        enum {
            KEYS_ONLY       = (Equals<Value, NullType>::VALUE),
            SCALE_FACTOR    = (CUB_MAX(sizeof(Key), sizeof(Value)) + 3) / 4,
            RADIX_BITS      = 5,
        };

        // UpsweepPolicy
        typedef BlockRadixSortUpsweepTilesPolicy <64,     CUB_MAX(1, 18 / SCALE_FACTOR), LOAD_LDG, RADIX_BITS> UpsweepPolicyKeys;
        typedef BlockRadixSortUpsweepTilesPolicy <128,    CUB_MAX(1, 15 / SCALE_FACTOR), LOAD_LDG, RADIX_BITS> UpsweepPolicyPairs;
        typedef typename If<KEYS_ONLY, UpsweepPolicyKeys, UpsweepPolicyPairs>::Type UpsweepPolicy;
/*
        // 4bit
        typedef BlockRadixSortUpsweepTilesPolicy <128, 15, LOAD_LDG, RADIX_BITS> UpsweepPolicyKeys;
        typedef BlockRadixSortUpsweepTilesPolicy <256, 13, LOAD_LDG, RADIX_BITS> UpsweepPolicyPairs;
*/
        // ScanPolicy
        typedef BlockScanTilesPolicy <1024, 4, BLOCK_LOAD_VECTORIZE, false, LOAD_DEFAULT, BLOCK_STORE_VECTORIZE, false, BLOCK_SCAN_RAKING_MEMOIZE> ScanPolicy;

        // DownsweepPolicy
        typedef BlockRadixSortDownsweepTilesPolicy <64,   CUB_MAX(1, 18 / SCALE_FACTOR), BLOCK_LOAD_DIRECT, LOAD_LDG, false, true, BLOCK_SCAN_WARP_SCANS, RADIX_SORT_SCATTER_TWO_PHASE, cudaSharedMemBankSizeEightByte, RADIX_BITS> DownsweepPolicyKeys;
        typedef BlockRadixSortDownsweepTilesPolicy <128,  CUB_MAX(1, 15 / SCALE_FACTOR), BLOCK_LOAD_DIRECT, LOAD_LDG, false, true, BLOCK_SCAN_WARP_SCANS, RADIX_SORT_SCATTER_TWO_PHASE, cudaSharedMemBankSizeEightByte, RADIX_BITS> DownsweepPolicyPairs;
        typedef typename If<KEYS_ONLY, DownsweepPolicyKeys, DownsweepPolicyPairs>::Type DownsweepPolicy;

/*
        // 4bit
        typedef BlockRadixSortDownsweepTilesPolicy <128, 15, BLOCK_LOAD_DIRECT, LOAD_LDG, false, true, BLOCK_SCAN_WARP_SCANS, RADIX_SORT_SCATTER_TWO_PHASE, cudaSharedMemBankSizeEightByte, RADIX_BITS> DownsweepPolicyKeys;
        typedef BlockRadixSortDownsweepTilesPolicy <256, 13, BLOCK_LOAD_DIRECT, LOAD_LDG, false, true, BLOCK_SCAN_WARP_SCANS, RADIX_SORT_SCATTER_TWO_PHASE, cudaSharedMemBankSizeEightByte, RADIX_BITS> DownsweepPolicyPairs;
*/
        enum { SUBSCRIPTION_FACTOR = 7 };
    };


    /// SM20 tune
    template <typename Key, typename Value, typename SizeT>
    struct TunedPolicies<Key, Value, SizeT, 200>
    {
        enum {
            KEYS_ONLY       = (Equals<Value, NullType>::VALUE),
            SCALE_FACTOR    = (CUB_MAX(sizeof(Key), sizeof(Value)) + 3) / 4,
            RADIX_BITS      = 5,
        };

        // UpsweepPolicy
        typedef BlockRadixSortUpsweepTilesPolicy <64, CUB_MAX(1, 18 / SCALE_FACTOR), LOAD_DEFAULT, RADIX_BITS> UpsweepPolicyKeys;
        typedef BlockRadixSortUpsweepTilesPolicy <128, CUB_MAX(1, 13 / SCALE_FACTOR), LOAD_DEFAULT, RADIX_BITS> UpsweepPolicyPairs;
        typedef typename If<KEYS_ONLY, UpsweepPolicyKeys, UpsweepPolicyPairs>::Type UpsweepPolicy;

        // ScanPolicy
        typedef BlockScanTilesPolicy <512, 4, BLOCK_LOAD_VECTORIZE, false, LOAD_DEFAULT, BLOCK_STORE_VECTORIZE, false, BLOCK_SCAN_RAKING_MEMOIZE> ScanPolicy;

        // DownsweepPolicy
        typedef BlockRadixSortDownsweepTilesPolicy <64, CUB_MAX(1, 18 / SCALE_FACTOR), BLOCK_LOAD_WARP_TRANSPOSE, LOAD_DEFAULT, false, false, BLOCK_SCAN_WARP_SCANS, RADIX_SORT_SCATTER_TWO_PHASE, cudaSharedMemBankSizeFourByte, RADIX_BITS> DownsweepPolicyKeys;
        typedef BlockRadixSortDownsweepTilesPolicy <128, CUB_MAX(1, 13 / SCALE_FACTOR), BLOCK_LOAD_WARP_TRANSPOSE, LOAD_DEFAULT, false, false, BLOCK_SCAN_WARP_SCANS, RADIX_SORT_SCATTER_TWO_PHASE, cudaSharedMemBankSizeFourByte, RADIX_BITS> DownsweepPolicyPairs;
        typedef typename If<KEYS_ONLY, DownsweepPolicyKeys, DownsweepPolicyPairs>::Type DownsweepPolicy;

        enum { SUBSCRIPTION_FACTOR = 3 };
    };


    /// SM10 tune
    template <typename Key, typename Value, typename SizeT>
    struct TunedPolicies<Key, Value, SizeT, 100>
    {
        enum {
            RADIX_BITS = 4,
        };

        // UpsweepPolicy
        typedef BlockRadixSortUpsweepTilesPolicy <64, 9, LOAD_DEFAULT, RADIX_BITS> UpsweepPolicy;

        // ScanPolicy
        typedef BlockScanTilesPolicy <256, 4, BLOCK_LOAD_VECTORIZE, false, LOAD_DEFAULT, BLOCK_STORE_VECTORIZE, false, BLOCK_SCAN_RAKING_MEMOIZE> ScanPolicy;

        // DownsweepPolicy
        typedef BlockRadixSortDownsweepTilesPolicy <64, 9, BLOCK_LOAD_WARP_TRANSPOSE, LOAD_DEFAULT, false, false, BLOCK_SCAN_WARP_SCANS, RADIX_SORT_SCATTER_TWO_PHASE, cudaSharedMemBankSizeFourByte, RADIX_BITS> DownsweepPolicy;

        enum { SUBSCRIPTION_FACTOR = 3 };
    };



    /******************************************************************************
     * Default policy initializer
     ******************************************************************************/

    /// Tuning policy for the PTX architecture that DeviceRadixSort operations will get dispatched to
    template <typename Key, typename Value, typename SizeT>
    struct PtxDefaultPolicies
    {

        static const int PTX_TUNE_ARCH =   (CUB_PTX_ARCH >= 350) ?
                                                350 :
                                                (CUB_PTX_ARCH >= 200) ?
                                                    200 :
                                                    100;

        // Tuned policy set for the current PTX compiler pass
        typedef TunedPolicies<Key, Value, SizeT, PTX_TUNE_ARCH> PtxTunedPolicies;

        // UpsweepPolicy that opaquely derives from the specialization corresponding to the current PTX compiler pass
        struct UpsweepPolicy : PtxTunedPolicies::UpsweepPolicy {};

        // ScanPolicy that opaquely derives from the specialization corresponding to the current PTX compiler pass
        struct ScanPolicy : PtxTunedPolicies::ScanPolicy {};

        // DownsweepPolicy that opaquely derives from the specialization corresponding to the current PTX compiler pass
        struct DownsweepPolicy : PtxTunedPolicies::DownsweepPolicy {};

        // Subscription factor for the current PTX compiler pass
        enum { SUBSCRIPTION_FACTOR = PtxTunedPolicies::SUBSCRIPTION_FACTOR };


        /**
         * Initialize dispatch params with the policies corresponding to the PTX assembly we will use
         */
        static void InitDispatchParams(
            int                    ptx_version,
            KernelDispachParams    &upsweep_dispatch_params,
            KernelDispachParams    &scan_dispatch_params,
            KernelDispachParams    &downsweep_dispatch_params)
        {
            if (ptx_version >= 350)
            {
                typedef TunedPolicies<Key, Value, SizeT, 350> TunedPolicies;
                upsweep_dispatch_params.InitUpsweepPolicy<typename TunedPolicies::UpsweepPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
                scan_dispatch_params.InitScanPolicy<typename TunedPolicies::ScanPolicy>();
                downsweep_dispatch_params.InitDownsweepPolicy<typename TunedPolicies::DownsweepPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
            }
            else if (ptx_version >= 200)
            {
                typedef TunedPolicies<Key, Value, SizeT, 200> TunedPolicies;
                upsweep_dispatch_params.InitUpsweepPolicy<typename TunedPolicies::UpsweepPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
                scan_dispatch_params.InitScanPolicy<typename TunedPolicies::ScanPolicy>();
                downsweep_dispatch_params.InitDownsweepPolicy<typename TunedPolicies::DownsweepPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
            }
            else
            {
                typedef TunedPolicies<Key, Value, SizeT, 100> TunedPolicies;
                upsweep_dispatch_params.InitUpsweepPolicy<typename TunedPolicies::UpsweepPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
                scan_dispatch_params.InitScanPolicy<typename TunedPolicies::ScanPolicy>();
                downsweep_dispatch_params.InitDownsweepPolicy<typename TunedPolicies::DownsweepPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
            }
        }
    };



    /******************************************************************************
     * Utility methods
     ******************************************************************************/

    /**
     * Internal dispatch routine for computing a device-wide reduction using a two-stages of kernel invocations.
     */
    template <
        typename            UpsweepKernelPtr,                       ///< Function type of cub::RadixSortUpsweepKernel
        typename            SpineKernelPtr,                         ///< Function type of cub::SpineScanKernel
        typename            DownsweepKernelPtr,                     ///< Function type of cub::RadixSortUpsweepKernel
        typename            Key,                                    ///< Key type
        typename            Value,                                  ///< Value type
        typename            SizeT>                                  ///< Integer type used for global array indexing
    __host__ __device__ __forceinline__
    static cudaError_t Dispatch(
        void                *d_temp_storage,                        ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                    ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        UpsweepKernelPtr    upsweep_kernel,                         ///< [in] Kernel function pointer to parameterization of cub::RadixSortUpsweepKernel
        SpineKernelPtr      scan_kernel,                            ///< [in] Kernel function pointer to parameterization of cub::SpineScanKernel
        DownsweepKernelPtr  downsweep_kernel,                       ///< [in] Kernel function pointer to parameterization of cub::RadixSortUpsweepKernel
        KernelDispachParams &upsweep_dispatch_params,               ///< [in] Dispatch parameters that match the policy that \p upsweep_kernel was compiled for
        KernelDispachParams &scan_dispatch_params,                  ///< [in] Dispatch parameters that match the policy that \p scan_kernel was compiled for
        KernelDispachParams &downsweep_dispatch_params,             ///< [in] Dispatch parameters that match the policy that \p downsweep_kernel was compiled for
        DoubleBuffer<Key>   &d_keys,                                ///< [in,out] Double-buffer whose current buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
        DoubleBuffer<Value> &d_values,                              ///< [in,out] Double-buffer whose current buffer contains the unsorted input values and, upon return, is updated to point to the sorted output values
        SizeT               num_items,                              ///< [in] Number of items to reduce
        int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The beginning (least-significant) bit index needed for key comparison
        int                 end_bit             = sizeof(Key) * 8,  ///< [in] <b>[optional]</b> The past-the-end (most-significant) bit index needed for key comparison
        cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
#ifndef CUB_RUNTIME_ENABLED

        // Kernel launch not supported from this device
        return CubDebug(cudaErrorNotSupported );

#else

        cudaError error = cudaSuccess;
        do
        {
            // Get device ordinal
            int device_ordinal;
            if (CubDebug(error = cudaGetDevice(&device_ordinal))) break;

            // Get SM count
            int sm_count;
            if (CubDebug(error = cudaDeviceGetAttribute (&sm_count, cudaDevAttrMultiProcessorCount, device_ordinal))) break;

            // Get a rough estimate of downsweep_kernel SM occupancy based upon the maximum SM occupancy of the targeted PTX architecture
            int downsweep_sm_occupancy = CUB_MIN(
                ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADBLOCKS,
                ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADS / downsweep_dispatch_params.block_threads);
            int upsweep_sm_occupancy = downsweep_sm_occupancy;

#ifndef __CUDA_ARCH__
            // We're on the host, so come up with more accurate estimates of SM occupancy from actual device properties
            Device device_props;
            if (CubDebug(error = device_props.Init(device_ordinal))) break;

            if (CubDebug(error = device_props.MaxSmOccupancy(
                downsweep_sm_occupancy,
                downsweep_kernel,
                downsweep_dispatch_params.block_threads))) break;

            if (CubDebug(error = device_props.MaxSmOccupancy(
                upsweep_sm_occupancy,
                upsweep_kernel,
                upsweep_dispatch_params.block_threads))) break;
#endif
            // Get device occupancies
            int downsweep_occupancy = downsweep_sm_occupancy * sm_count;

            // Get even-share work distribution descriptor
            GridEvenShare<SizeT> even_share;
            int max_downsweep_grid_size = downsweep_occupancy * downsweep_dispatch_params.subscription_factor;
            int downsweep_grid_size;
            even_share.GridInit(num_items, max_downsweep_grid_size, downsweep_dispatch_params.tile_size);
            downsweep_grid_size = even_share.grid_size;

            // Get number of spine elements (round up to nearest spine scan kernel tile size)
            int bins            = 1 << downsweep_dispatch_params.radix_bits;
            int spine_size      = downsweep_grid_size * bins;
            int spine_tiles     = (spine_size + scan_dispatch_params.tile_size - 1) / scan_dispatch_params.tile_size;
            spine_size          = spine_tiles * scan_dispatch_params.tile_size;

            int alt_bins            = 1 << downsweep_dispatch_params.alt_radix_bits;
            int alt_spine_size      = downsweep_grid_size * alt_bins;
            int alt_spine_tiles     = (alt_spine_size + scan_dispatch_params.tile_size - 1) / scan_dispatch_params.tile_size;
            alt_spine_size          = alt_spine_tiles * scan_dispatch_params.tile_size;

            // Temporary storage allocation requirements
            void* allocations[1];
            size_t allocation_sizes[1] =
            {
                spine_size * sizeof(SizeT),    // bytes needed for privatized block digit histograms
            };

            // Alias temporaries (or set the necessary size of the storage allocation)
            if (CubDebug(error = AliasTemporaries(d_temp_storage, temp_storage_bytes, allocations, allocation_sizes))) break;

            // Return if the caller is simply requesting the size of the storage allocation
            if (d_temp_storage == NULL)
                return cudaSuccess;

            // Privatized per-block digit histograms
            SizeT *d_spine = (SizeT*) allocations[0];

#ifndef __CUDA_ARCH__
            // Get current smem bank configuration
            cudaSharedMemConfig original_smem_config;
            if (CubDebug(error = cudaDeviceGetSharedMemConfig(&original_smem_config))) break;
            cudaSharedMemConfig current_smem_config = original_smem_config;
#endif
            // Iterate over digit places
            int current_bit = begin_bit;
            while (current_bit < end_bit)
            {
                // Use primary bit granularity if bits remaining is a whole multiple of bit primary granularity
                int bits_remaining = end_bit - current_bit;
                bool use_primary_bit_granularity = (bits_remaining % downsweep_dispatch_params.radix_bits == 0);
                int radix_bits = (use_primary_bit_granularity) ?
                    downsweep_dispatch_params.radix_bits :
                    downsweep_dispatch_params.alt_radix_bits;

#ifndef __CUDA_ARCH__
                // Update smem config if necessary
                if (current_smem_config != upsweep_dispatch_params.smem_config)
                {
                    if (CubDebug(error = cudaDeviceSetSharedMemConfig(upsweep_dispatch_params.smem_config))) break;
                    current_smem_config = upsweep_dispatch_params.smem_config;
                }
#endif

                // Log upsweep_kernel configuration
                if (stream_synchronous)
                    CubLog("Invoking upsweep_kernel<<<%d, %d, 0, %lld>>>(), %d smem config, %d items per thread, %d SM occupancy, selector %d, current bit %d, bit_grain %d\n",
                    downsweep_grid_size, upsweep_dispatch_params.block_threads, (long long) stream, upsweep_dispatch_params.smem_config, upsweep_dispatch_params.items_per_thread, upsweep_sm_occupancy, d_keys.selector, current_bit, radix_bits);

                // Invoke upsweep_kernel with same grid size as downsweep_kernel
                upsweep_kernel<<<downsweep_grid_size, upsweep_dispatch_params.block_threads, 0, stream>>>(
                    d_keys.d_buffers[d_keys.selector],
                    d_spine,
                    num_items,
                    current_bit,
                    use_primary_bit_granularity,
                    (current_bit == begin_bit),
                    even_share);

                // Sync the stream if specified
                if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;

                // Log scan_kernel configuration
                if (stream_synchronous) CubLog("Invoking scan_kernel<<<%d, %d, 0, %lld>>>(), %d items per thread\n",
                    1, scan_dispatch_params.block_threads, (long long) stream, scan_dispatch_params.items_per_thread);

                // Invoke scan_kernel
                scan_kernel<<<1, scan_dispatch_params.block_threads, 0, stream>>>(
                    d_spine,
                    (use_primary_bit_granularity) ? spine_size : alt_spine_size);

                // Sync the stream if specified
                if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;

#ifndef __CUDA_ARCH__
                // Update smem config if necessary
                if (current_smem_config != downsweep_dispatch_params.smem_config)
                {
                    if (CubDebug(error = cudaDeviceSetSharedMemConfig(downsweep_dispatch_params.smem_config))) break;
                    current_smem_config = downsweep_dispatch_params.smem_config;
                }
#endif

                // Log downsweep_kernel configuration
                if (stream_synchronous) CubLog("Invoking downsweep_kernel<<<%d, %d, 0, %lld>>>(), %d smem config, %d items per thread, %d SM occupancy\n",
                    downsweep_grid_size, downsweep_dispatch_params.block_threads, (long long) stream, downsweep_dispatch_params.smem_config, downsweep_dispatch_params.items_per_thread, downsweep_sm_occupancy);

                // Invoke downsweep_kernel
                downsweep_kernel<<<downsweep_grid_size, downsweep_dispatch_params.block_threads, 0, stream>>>(
                    d_keys.d_buffers[d_keys.selector],
                    d_keys.d_buffers[d_keys.selector ^ 1],
                    d_values.d_buffers[d_values.selector],
                    d_values.d_buffers[d_values.selector ^ 1],
                    d_spine,
                    num_items,
                    current_bit,
                    use_primary_bit_granularity,
                    (current_bit == begin_bit),
                    (current_bit + downsweep_dispatch_params.radix_bits >= end_bit),
                    even_share);

                // Sync the stream if specified
                if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;

                // Invert selectors
                d_keys.selector ^= 1;
                d_values.selector ^= 1;

                // Update current bit position
                current_bit += radix_bits;
            }

#ifndef __CUDA_ARCH__
            // Reset smem config if necessary
            if (current_smem_config != original_smem_config)
            {
                if (CubDebug(error = cudaDeviceSetSharedMemConfig(original_smem_config))) break;
            }
#endif

        }
        while (0);

        return error;

#endif // CUB_RUNTIME_ENABLED
    }



    #endif // DOXYGEN_SHOULD_SKIP_THIS

    /******************************************************************************
     * Interface
     ******************************************************************************/


    /**
     * \brief Sorts key-value pairs.
     *
     * \par
     * The sorting operation requires a pair of key buffers and a pair of value
     * buffers.  Each pair is wrapped in a DoubleBuffer structure whose member
     * DoubleBuffer::Current() references the active buffer.  The currently-active
     * buffer may be changed by the sorting operation.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \par
     * The code snippet below illustrates the sorting of a device vector of \p int keys
     * with associated vector of \p int values.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Create a set of DoubleBuffers to wrap pairs of device pointers for
     * // sorting data (keys, values, and equivalently-sized alternate buffers)
     * int num_items = ...
     * cub::DoubleBuffer<int> d_keys(d_key_buf, d_key_alt_buf);
     * cub::DoubleBuffer<int> d_values(d_value_buf, d_value_alt_buf);
     *
     * // Determine temporary device storage requirements for sorting operation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, num_items);
     *
     * // Allocate temporary storage for sorting operation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Run sorting operation
     * cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, num_items);
     *
     * // Sorted keys and values are referenced by d_keys.Current() and d_values.Current()
     *
     * \endcode
     *
     * \tparam Key      <b>[inferred]</b> Key type
     * \tparam Value    <b>[inferred]</b> Value type
     */
    template <
        typename            Key,
        typename            Value>
    __host__ __device__ __forceinline__
    static cudaError_t SortPairs(
        void                *d_temp_storage,                        ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                    ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        DoubleBuffer<Key>   &d_keys,                                ///< [in,out] Double-buffer of keys whose current buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
        DoubleBuffer<Value> &d_values,                              ///< [in,out] Double-buffer of values whose current buffer contains the unsorted input values and, upon return, is updated to point to the sorted output values
        int                 num_items,                              ///< [in] Number of items to reduce
        int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The first (least-significant) bit index needed for key comparison
        int                 end_bit             = sizeof(Key) * 8,  ///< [in] <b>[optional]</b> The past-the-end (most-significant) bit index needed for key comparison
        cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
        // Type used for array indexing
        typedef int SizeT;

        // Tuning polices
        typedef PtxDefaultPolicies<Key, Value, SizeT>           PtxDefaultPolicies; // Wrapper of default kernel policies
        typedef typename PtxDefaultPolicies::UpsweepPolicy      UpsweepPolicy;      // Upsweep kernel policy
        typedef typename PtxDefaultPolicies::ScanPolicy         ScanPolicy;         // Scan kernel policy
        typedef typename PtxDefaultPolicies::DownsweepPolicy    DownsweepPolicy;    // Downsweep kernel policy

        cudaError error = cudaSuccess;
        do
        {
            // Declare dispatch parameters
            KernelDispachParams upsweep_dispatch_params;
            KernelDispachParams scan_dispatch_params;
            KernelDispachParams downsweep_dispatch_params;

#ifdef __CUDA_ARCH__
            // We're on the device, so initialize the dispatch parameters with the PtxDefaultPolicies directly
            upsweep_dispatch_params.InitUpsweepPolicy<UpsweepPolicy>(PtxDefaultPolicies::SUBSCRIPTION_FACTOR);
            scan_dispatch_params.InitScanPolicy<ScanPolicy>();
            downsweep_dispatch_params.InitDownsweepPolicy<DownsweepPolicy>(PtxDefaultPolicies::SUBSCRIPTION_FACTOR);
#else
            // We're on the host, so lookup and initialize the dispatch parameters with the policies that match the device's PTX version
            int ptx_version;
            if (CubDebug(error = PtxVersion(ptx_version))) break;
            PtxDefaultPolicies::InitDispatchParams(
                ptx_version,
                upsweep_dispatch_params,
                scan_dispatch_params,
                downsweep_dispatch_params);
#endif
            // Dispatch
            if (CubDebug(error = Dispatch(
                d_temp_storage,
                temp_storage_bytes,
                RadixSortUpsweepKernel<UpsweepPolicy, Key, SizeT>,
                RadixSortScanKernel<ScanPolicy, SizeT>,
                RadixSortDownsweepKernel<DownsweepPolicy, Key, Value, SizeT>,
                upsweep_dispatch_params,
                scan_dispatch_params,
                downsweep_dispatch_params,
                d_keys,
                d_values,
                num_items,
                begin_bit,
                end_bit,
                stream,
                stream_synchronous))) break;
        }
        while (0);

        return error;
    }


    /**
     * \brief Sorts keys
     *
     * \par
     * The sorting operation requires a pair of key buffers.  The pair is
     * wrapped in a DoubleBuffer structure whose member DoubleBuffer::Current()
     * references the active buffer.  The currently-active buffer may be changed
     * by the sorting operation.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \par
     * The code snippet below illustrates the sorting of a device vector of \p int keys.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Create a set of DoubleBuffers to wrap pairs of device pointers for
     * // sorting data (keys and equivalently-sized alternate buffer)
     * int num_items = ...
     * cub::DoubleBuffer<int> d_keys(d_key_buf, d_key_alt_buf);
     *
     * // Determine temporary device storage requirements for sorting operation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, num_items);
     *
     * // Allocate temporary storage for sorting operation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Run sorting operation
     * cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, num_items);
     *
     * // Sorted keys are referenced by d_keys.Current()
     *
     * \endcode
     *
     * \tparam Key      <b>[inferred]</b> Key type
     */
    template <typename Key>
    __host__ __device__ __forceinline__
    static cudaError_t SortKeys(
        void                *d_temp_storage,                        ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                    ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        DoubleBuffer<Key>   &d_keys,                                ///< [in,out] Double-buffer of keys whose current buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
        int                 num_items,                              ///< [in] Number of items to reduce
        int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The first (least-significant) bit index needed for key comparison
        int                 end_bit             = sizeof(Key) * 8,  ///< [in] <b>[optional]</b> The past-the-end (most-significant) bit index needed for key comparison
        cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
        DoubleBuffer<NullType> d_values;
        return SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, num_items, begin_bit, end_bit, stream, stream_synchronous);
    }

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)


