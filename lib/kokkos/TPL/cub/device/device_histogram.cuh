
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
 * cub::DeviceHistogram provides device-wide parallel operations for constructing histogram(s) from samples data residing within global memory.
 */

#pragma once

#include <stdio.h>
#include <iterator>

#include "block/block_histo_tiles.cuh"
#include "../grid/grid_even_share.cuh"
#include "../grid/grid_queue.cuh"
#include "../util_debug.cuh"
#include "../util_device.cuh"
#include "../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/******************************************************************************
 * Kernel entry points
 *****************************************************************************/

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document


/**
 * Initialization pass kernel entry point (multi-block).  Prepares queue descriptors zeroes global counters.
 */
template <
    int                                             BINS,                   ///< Number of histogram bins per channel
    int                                             ACTIVE_CHANNELS,        ///< Number of channels actively being histogrammed
    typename                                        SizeT,                  ///< Integer type used for global array indexing
    typename                                        HistoCounter>           ///< Integral type for counting sample occurrences per histogram bin
__launch_bounds__ (BINS, 1)
__global__ void InitHistoKernel(
    GridQueue<SizeT>                                grid_queue,             ///< [in] Descriptor for performing dynamic mapping of tile data to thread blocks
    ArrayWrapper<HistoCounter*, ACTIVE_CHANNELS>    d_out_histograms,       ///< [out] Histogram counter data having logical dimensions <tt>HistoCounter[ACTIVE_CHANNELS][BINS]</tt>
    SizeT                                           num_samples)            ///< [in] Total number of samples \p d_samples for all channels
{
    d_out_histograms.array[blockIdx.x][threadIdx.x] = 0;
    if (threadIdx.x == 0) grid_queue.ResetDrain(num_samples);
}


/**
 * Histogram pass kernel entry point (multi-block).  Computes privatized histograms, one per thread block.
 */
template <
    typename                                        BlockHistogramTilesPolicy,   ///< Tuning policy for cub::BlockHistogramTiles abstraction
    int                                             BINS,                       ///< Number of histogram bins per channel
    int                                             CHANNELS,                   ///< Number of channels interleaved in the input data (may be greater than the number of channels being actively histogrammed)
    int                                             ACTIVE_CHANNELS,            ///< Number of channels actively being histogrammed
    typename                                        InputIteratorRA,            ///< The input iterator type (may be a simple pointer type).  Must have a value type that is assignable to <tt>unsigned char</tt>
    typename                                        HistoCounter,               ///< Integral type for counting sample occurrences per histogram bin
    typename                                        SizeT>                      ///< Integer type used for global array indexing
__launch_bounds__ (int(BlockHistogramTilesPolicy::BLOCK_THREADS), BlockHistogramTilesPolicy::SM_OCCUPANCY)
__global__ void MultiBlockHistogramKernel(
    InputIteratorRA                                 d_samples,                  ///< [in] Array of sample data. The samples from different channels are assumed to be interleaved (e.g., an array of 32b pixels where each pixel consists of four RGBA 8b samples).
    ArrayWrapper<HistoCounter*, ACTIVE_CHANNELS>    d_out_histograms,           ///< [out] Histogram counter data having logical dimensions <tt>HistoCounter[ACTIVE_CHANNELS][gridDim.x][BINS]</tt>
    SizeT                                           num_samples,                ///< [in] Total number of samples \p d_samples for all channels
    GridEvenShare<SizeT>                            even_share,                 ///< [in] Descriptor for how to map an even-share of tiles across thread blocks
    GridQueue<SizeT>                                queue)                      ///< [in] Descriptor for performing dynamic mapping of tile data to thread blocks
{
    // Constants
    enum
    {
        BLOCK_THREADS       = BlockHistogramTilesPolicy::BLOCK_THREADS,
        ITEMS_PER_THREAD    = BlockHistogramTilesPolicy::ITEMS_PER_THREAD,
        TILE_SIZE           = BLOCK_THREADS * ITEMS_PER_THREAD,
    };

    // Thread block type for compositing input tiles
    typedef BlockHistogramTiles<BlockHistogramTilesPolicy, BINS, CHANNELS, ACTIVE_CHANNELS, InputIteratorRA, HistoCounter, SizeT> BlockHistogramTilesT;

    // Shared memory for BlockHistogramTiles
    __shared__ typename BlockHistogramTilesT::TempStorage temp_storage;

    // Consume input tiles
    BlockHistogramTilesT(temp_storage, d_samples, d_out_histograms.array).ConsumeTiles(
        num_samples,
        even_share,
        queue,
        Int2Type<BlockHistogramTilesPolicy::GRID_MAPPING>());
}


/**
 * Block-aggregation pass kernel entry point (single-block).  Aggregates privatized threadblock histograms from a previous multi-block histogram pass.
 */
template <
    int                                             BINS,                   ///< Number of histogram bins per channel
    int                                             ACTIVE_CHANNELS,        ///< Number of channels actively being histogrammed
    typename                                        HistoCounter>           ///< Integral type for counting sample occurrences per histogram bin
__launch_bounds__ (BINS, 1)
__global__ void AggregateHistoKernel(
    HistoCounter*                                   d_block_histograms,     ///< [in] Histogram counter data having logical dimensions <tt>HistoCounter[ACTIVE_CHANNELS][num_threadblocks][BINS]</tt>
    ArrayWrapper<HistoCounter*, ACTIVE_CHANNELS>    d_out_histograms,       ///< [out] Histogram counter data having logical dimensions <tt>HistoCounter[ACTIVE_CHANNELS][BINS]</tt>
    int                                             num_threadblocks)       ///< [in] Number of threadblock histograms per channel in \p d_block_histograms
{
    // Accumulate threadblock-histograms from the channel
    HistoCounter bin_aggregate = 0;

    int block_offset = blockIdx.x * (num_threadblocks * BINS);
    int block_oob = block_offset + (num_threadblocks * BINS);

#if CUB_PTX_ARCH >= 200
    #pragma unroll 32
#endif
    while (block_offset < block_oob)
    {
        bin_aggregate += d_block_histograms[block_offset + threadIdx.x];
        block_offset += BINS;
    }

    // Output
    d_out_histograms.array[blockIdx.x][threadIdx.x] = bin_aggregate;
}

#endif // DOXYGEN_SHOULD_SKIP_THIS



/******************************************************************************
 * DeviceHistogram
 *****************************************************************************/

/**
 * \brief DeviceHistogram provides device-wide parallel operations for constructing histogram(s) from samples data residing within global memory. ![](histogram_logo.png)
 * \ingroup DeviceModule
 *
 * \par Overview
 * A <a href="http://en.wikipedia.org/wiki/Histogram"><em>histogram</em></a>
 * counts the number of observations that fall into each of the disjoint categories (known as <em>bins</em>).
 *
 * \par Usage Considerations
 * \cdp_class{DeviceHistogram}
 *
 * \par Performance
 *
 * \image html histo_perf.png
 *
 */
struct DeviceHistogram
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    /// Generic structure for encapsulating dispatch properties.  Mirrors the constants within BlockHistogramTilesPolicy.
    struct KernelDispachParams
    {
        // Policy fields
        int                         block_threads;
        int                         items_per_thread;
        BlockHistogramTilesAlgorithm    block_algorithm;
        GridMappingStrategy         grid_mapping;
        int                         subscription_factor;

        // Derived fields
        int                         channel_tile_size;

        template <typename BlockHistogramTilesPolicy>
        __host__ __device__ __forceinline__
        void Init(int subscription_factor = 1)
        {
            block_threads               = BlockHistogramTilesPolicy::BLOCK_THREADS;
            items_per_thread            = BlockHistogramTilesPolicy::ITEMS_PER_THREAD;
            block_algorithm             = BlockHistogramTilesPolicy::GRID_ALGORITHM;
            grid_mapping                = BlockHistogramTilesPolicy::GRID_MAPPING;
            this->subscription_factor   = subscription_factor;

            channel_tile_size           = block_threads * items_per_thread;
        }

        __host__ __device__ __forceinline__
        void Print()
        {
            printf("%d, %d, %d, %d, %d",
                block_threads,
                items_per_thread,
                block_algorithm,
                grid_mapping,
                subscription_factor);
        }

    };


    /******************************************************************************
     * Tuning policies
     ******************************************************************************/

    /// Specializations of tuned policy types for different PTX architectures
    template <
        int                         CHANNELS,
        int                         ACTIVE_CHANNELS,
        BlockHistogramTilesAlgorithm    GRID_ALGORITHM,
        int                         ARCH>
    struct TunedPolicies;

    /// SM35 tune
    template <int CHANNELS, int ACTIVE_CHANNELS, BlockHistogramTilesAlgorithm GRID_ALGORITHM>
    struct TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 350>
    {
        typedef BlockHistogramTilesPolicy<
            (GRID_ALGORITHM == GRID_HISTO_SORT) ? 128 : 256,
            (GRID_ALGORITHM == GRID_HISTO_SORT) ? 12 : (30 / ACTIVE_CHANNELS),
            GRID_ALGORITHM,
            (GRID_ALGORITHM == GRID_HISTO_SORT) ? GRID_MAPPING_DYNAMIC : GRID_MAPPING_EVEN_SHARE,
            (GRID_ALGORITHM == GRID_HISTO_SORT) ? 8 : 1> MultiBlockPolicy;
        enum { SUBSCRIPTION_FACTOR = 7 };
    };

    /// SM30 tune
    template <int CHANNELS, int ACTIVE_CHANNELS, BlockHistogramTilesAlgorithm GRID_ALGORITHM>
    struct TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 300>
    {
        typedef BlockHistogramTilesPolicy<
            128,
            (GRID_ALGORITHM == GRID_HISTO_SORT) ? 20 : (22 / ACTIVE_CHANNELS),
            GRID_ALGORITHM,
            (GRID_ALGORITHM == GRID_HISTO_SORT) ? GRID_MAPPING_DYNAMIC : GRID_MAPPING_EVEN_SHARE,
            1> MultiBlockPolicy;
        enum { SUBSCRIPTION_FACTOR = 1 };
    };

    /// SM20 tune
    template <int CHANNELS, int ACTIVE_CHANNELS, BlockHistogramTilesAlgorithm GRID_ALGORITHM>
    struct TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 200>
    {
        typedef BlockHistogramTilesPolicy<
            128,
            (GRID_ALGORITHM == GRID_HISTO_SORT) ? 21 : (23 / ACTIVE_CHANNELS),
            GRID_ALGORITHM,
            GRID_MAPPING_DYNAMIC,
            1> MultiBlockPolicy;
        enum { SUBSCRIPTION_FACTOR = 1 };
    };

    /// SM10 tune
    template <int CHANNELS, int ACTIVE_CHANNELS, BlockHistogramTilesAlgorithm GRID_ALGORITHM>
    struct TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 100>
    {
        typedef BlockHistogramTilesPolicy<
            128, 
            7, 
            GRID_HISTO_SORT,        // (use sort regardless because atomics are perf-useless)
            GRID_MAPPING_EVEN_SHARE,
            1> MultiBlockPolicy;
        enum { SUBSCRIPTION_FACTOR = 1 };
    };


    /// Tuning policy for the PTX architecture that DeviceHistogram operations will get dispatched to
    template <
        int                         CHANNELS,
        int                         ACTIVE_CHANNELS,
        BlockHistogramTilesAlgorithm      GRID_ALGORITHM>
    struct PtxDefaultPolicies
    {
        static const int PTX_TUNE_ARCH =   (CUB_PTX_ARCH >= 350) ?
                                                350 :
                                                (CUB_PTX_ARCH >= 300) ?
                                                    300 :
                                                    (CUB_PTX_ARCH >= 200) ?
                                                        200 :
                                                        100;

        // Tuned policy set for the current PTX compiler pass
        typedef TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, PTX_TUNE_ARCH> PtxTunedPolicies;

        // Subscription factor for the current PTX compiler pass
        static const int SUBSCRIPTION_FACTOR = PtxTunedPolicies::SUBSCRIPTION_FACTOR;

        // MultiBlockPolicy that opaquely derives from the specialization corresponding to the current PTX compiler pass
        struct MultiBlockPolicy : PtxTunedPolicies::MultiBlockPolicy {};

        /**
         * Initialize dispatch params with the policies corresponding to the PTX assembly we will use
         */
        static void InitDispatchParams(int ptx_version, KernelDispachParams &multi_block_dispatch_params)
        {
            if (ptx_version >= 350)
            {
                typedef TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 350> TunedPolicies;
                multi_block_dispatch_params.Init<typename TunedPolicies::MultiBlockPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
            }
            else if (ptx_version >= 300)
            {
                typedef TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 300> TunedPolicies;
                multi_block_dispatch_params.Init<typename TunedPolicies::MultiBlockPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
            }
            else if (ptx_version >= 200)
            {
                typedef TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 200> TunedPolicies;
                multi_block_dispatch_params.Init<typename TunedPolicies::MultiBlockPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
            }
            else
            {
                typedef TunedPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM, 100> TunedPolicies;
                multi_block_dispatch_params.Init<typename TunedPolicies::MultiBlockPolicy>(TunedPolicies::SUBSCRIPTION_FACTOR);
            }
        }
    };


    /******************************************************************************
     * Utility methods
     ******************************************************************************/

    /**
     * Internal dispatch routine for invoking device-wide, multi-channel, histogram
     */
    template <
        int                         BINS,                               ///< Number of histogram bins per channel
        int                         CHANNELS,                           ///< Number of channels interleaved in the input data (may be greater than the number of channels being actively histogrammed)
        int                         ACTIVE_CHANNELS,                    ///< Number of channels actively being histogrammed
        typename                    InitHistoKernelPtr,                 ///< Function type of cub::InitHistoKernel
        typename                    MultiBlockHistogramKernelPtr,           ///< Function type of cub::MultiBlockHistogramKernel
        typename                    AggregateHistoKernelPtr,            ///< Function type of cub::AggregateHistoKernel
        typename                    InputIteratorRA,                    ///< The input iterator type (may be a simple pointer type).  Must have a value type that is assignable to <tt>unsigned char</tt>
        typename                    HistoCounter,                       ///< Integral type for counting sample occurrences per histogram bin
        typename                    SizeT>                              ///< Integer type used for global array indexing
    __host__ __device__ __forceinline__
    static cudaError_t Dispatch(
        void                        *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t                      &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InitHistoKernelPtr          init_kernel,                        ///< [in] Kernel function pointer to parameterization of cub::InitHistoKernel
        MultiBlockHistogramKernelPtr    multi_block_kernel,                 ///< [in] Kernel function pointer to parameterization of cub::MultiBlockHistogramKernel
        AggregateHistoKernelPtr     aggregate_kernel,                   ///< [in] Kernel function pointer to parameterization of cub::AggregateHistoKernel
        KernelDispachParams         &multi_block_dispatch_params,       ///< [in] Dispatch parameters that match the policy that \p multi_block_kernel was compiled for
        InputIteratorRA             d_samples,                          ///< [in] Input samples to histogram
        HistoCounter                *d_histograms[ACTIVE_CHANNELS],     ///< [out] Array of channel histograms, each having BINS counters of integral type \p HistoCounter.
        SizeT                       num_samples,                        ///< [in] Number of samples to process
        cudaStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                        stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
#ifndef CUB_RUNTIME_ENABLED

        // Kernel launch not supported from this device
        return CubDebug(cudaErrorNotSupported);

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

            // Get a rough estimate of multi_block_kernel SM occupancy based upon the maximum SM occupancy of the targeted PTX architecture
            int multi_block_sm_occupancy = CUB_MIN(
                ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADBLOCKS,
                ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADS / multi_block_dispatch_params.block_threads);

#ifndef __CUDA_ARCH__
            // We're on the host, so come up with a more accurate estimate of multi_block_kernel SM occupancy from actual device properties
            Device device_props;
            if (CubDebug(error = device_props.Init(device_ordinal))) break;

            if (CubDebug(error = device_props.MaxSmOccupancy(
                multi_block_sm_occupancy,
                multi_block_kernel,
                multi_block_dispatch_params.block_threads))) break;
#endif

            // Get device occupancy for multi_block_kernel
            int multi_block_occupancy = multi_block_sm_occupancy * sm_count;

            // Even-share work distribution
            GridEvenShare<SizeT> even_share;

            // Get tile size for multi_block_kernel
            int multi_block_tile_size = multi_block_dispatch_params.channel_tile_size * CHANNELS;

            // Get grid size for multi_block_kernel
            int multi_block_grid_size;
            switch (multi_block_dispatch_params.grid_mapping)
            {
            case GRID_MAPPING_EVEN_SHARE:

                // Work is distributed evenly
                even_share.GridInit(
                    num_samples,
                    multi_block_occupancy * multi_block_dispatch_params.subscription_factor,
                    multi_block_tile_size);
                multi_block_grid_size = even_share.grid_size;
                break;

            case GRID_MAPPING_DYNAMIC:

                // Work is distributed dynamically
                int num_tiles           = (num_samples + multi_block_tile_size - 1) / multi_block_tile_size;
                multi_block_grid_size   = (num_tiles < multi_block_occupancy) ?
                    num_tiles :                 // Not enough to fill the device with threadblocks
                    multi_block_occupancy;      // Fill the device with threadblocks
                break;
            };

            // Temporary storage allocation requirements
            void* allocations[2];
            size_t allocation_sizes[2] =
            {
                ACTIVE_CHANNELS * multi_block_grid_size * sizeof(HistoCounter) * BINS,      // bytes needed for privatized histograms
                GridQueue<int>::AllocationSize()                                            // bytes needed for grid queue descriptor
            };

            if (CubDebug(error = AliasTemporaries(d_temp_storage, temp_storage_bytes, allocations, allocation_sizes))) break;

            // Return if the caller is simply requesting the size of the storage allocation
            if (d_temp_storage == NULL)
                return cudaSuccess;

            // Privatized per-block reductions
            HistoCounter *d_block_histograms = (HistoCounter*) allocations[0];

            // Grid queue descriptor
            GridQueue<SizeT> queue(allocations[1]);

            // Setup array wrapper for histogram channel output (because we can't pass static arrays as kernel parameters)
            ArrayWrapper<HistoCounter*, ACTIVE_CHANNELS> d_histo_wrapper;
            for (int CHANNEL = 0; CHANNEL < ACTIVE_CHANNELS; ++CHANNEL)
                d_histo_wrapper.array[CHANNEL] = d_histograms[CHANNEL];

            // Setup array wrapper for temporary histogram channel output (because we can't pass static arrays as kernel parameters)
            ArrayWrapper<HistoCounter*, ACTIVE_CHANNELS> d_temp_histo_wrapper;
            for (int CHANNEL = 0; CHANNEL < ACTIVE_CHANNELS; ++CHANNEL)
                d_temp_histo_wrapper.array[CHANNEL] = d_block_histograms + (CHANNEL * multi_block_grid_size * BINS);

            // Log init_kernel configuration
            if (stream_synchronous) CubLog("Invoking init_kernel<<<%d, %d, 0, %lld>>>()\n", ACTIVE_CHANNELS, BINS, (long long) stream);

            // Invoke init_kernel to initialize counters and queue descriptor
            init_kernel<<<ACTIVE_CHANNELS, BINS, 0, stream>>>(queue, d_histo_wrapper, num_samples);

            // Sync the stream if specified
            if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;

            // Whether we need privatized histograms (i.e., non-global atomics and multi-block)
            bool privatized_temporaries = (multi_block_grid_size > 1) && (multi_block_dispatch_params.block_algorithm != GRID_HISTO_GLOBAL_ATOMIC);

            // Log multi_block_kernel configuration
            if (stream_synchronous) CubLog("Invoking multi_block_kernel<<<%d, %d, 0, %lld>>>(), %d items per thread, %d SM occupancy\n",
                multi_block_grid_size, multi_block_dispatch_params.block_threads, (long long) stream, multi_block_dispatch_params.items_per_thread, multi_block_sm_occupancy);

            // Invoke multi_block_kernel
            multi_block_kernel<<<multi_block_grid_size, multi_block_dispatch_params.block_threads, 0, stream>>>(
                d_samples,
                (privatized_temporaries) ?
                    d_temp_histo_wrapper :
                    d_histo_wrapper,
                num_samples,
                even_share,
                queue);

            // Sync the stream if specified
            if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;

            // Aggregate privatized block histograms if necessary
            if (privatized_temporaries)
            {
                // Log aggregate_kernel configuration
                if (stream_synchronous) CubLog("Invoking aggregate_kernel<<<%d, %d, 0, %lld>>>()\n",
                    ACTIVE_CHANNELS, BINS, (long long) stream);

                // Invoke aggregate_kernel
                aggregate_kernel<<<ACTIVE_CHANNELS, BINS, 0, stream>>>(
                    d_block_histograms,
                    d_histo_wrapper,
                    multi_block_grid_size);

                // Sync the stream if specified
                if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;
            }
        }
        while (0);

        return error;
#endif // CUB_RUNTIME_ENABLED
    }


    /**
     * \brief Computes a device-wide histogram
     *
     * \tparam GRID_ALGORITHM      cub::BlockHistogramTilesAlgorithm enumerator specifying the underlying algorithm to use
     * \tparam CHANNELS             Number of channels interleaved in the input data (may be greater than the number of channels being actively histogrammed)
     * \tparam ACTIVE_CHANNELS      <b>[inferred]</b> Number of channels actively being histogrammed
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)  Must have a value type that is assignable to <tt>unsigned char</tt>
     * \tparam HistoCounter         <b>[inferred]</b> Integral type for counting sample occurrences per histogram bin
     */
    template <
        BlockHistogramTilesAlgorithm    GRID_ALGORITHM,
        int                         BINS,                       ///< Number of histogram bins per channel
        int                         CHANNELS,                   ///< Number of channels interleaved in the input data (may be greater than the number of channels being actively histogrammed)
        int                         ACTIVE_CHANNELS,            ///< Number of channels actively being histogrammed
        typename                    InputIteratorRA,            ///< The input iterator type (may be a simple pointer type).  Must have a value type that is assignable to <tt>unsigned char</tt>
        typename                    HistoCounter>               ///< Integral type for counting sample occurrences per histogram bin
    __host__ __device__ __forceinline__
    static cudaError_t Dispatch(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_samples,                          ///< [in] Input samples to histogram
        HistoCounter        *d_histograms[ACTIVE_CHANNELS],     ///< [out] Array of channel histograms, each having BINS counters of integral type \p HistoCounter.
        int                 num_samples,                        ///< [in] Number of samples to process
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
        // Type used for array indexing
        typedef int SizeT;

        // Tuning polices for the PTX architecture that will get dispatched to
        typedef PtxDefaultPolicies<CHANNELS, ACTIVE_CHANNELS, GRID_ALGORITHM> PtxDefaultPolicies;
        typedef typename PtxDefaultPolicies::MultiBlockPolicy MultiBlockPolicy;

        cudaError error = cudaSuccess;
        do
        {
            // Declare dispatch parameters
            KernelDispachParams multi_block_dispatch_params;

        #ifdef __CUDA_ARCH__

            // We're on the device, so initialize the dispatch parameters with the PtxDefaultPolicies directly
            multi_block_dispatch_params.Init<MultiBlockPolicy>(PtxDefaultPolicies::SUBSCRIPTION_FACTOR);

        #else

            // We're on the host, so lookup and initialize the dispatch parameters with the policies that match the device's PTX version
            int ptx_version;
            if (CubDebug(error = PtxVersion(ptx_version))) break;
            PtxDefaultPolicies::InitDispatchParams(ptx_version, multi_block_dispatch_params);

        #endif

            Dispatch<BINS, CHANNELS, ACTIVE_CHANNELS>(
                d_temp_storage,
                temp_storage_bytes,
                InitHistoKernel<BINS, ACTIVE_CHANNELS, SizeT, HistoCounter>,
                MultiBlockHistogramKernel<MultiBlockPolicy, BINS, CHANNELS, ACTIVE_CHANNELS, InputIteratorRA, HistoCounter, SizeT>,
                AggregateHistoKernel<BINS, ACTIVE_CHANNELS, HistoCounter>,
                multi_block_dispatch_params,
                d_samples,
                d_histograms,
                num_samples,
                stream,
                stream_synchronous);

            if (CubDebug(error)) break;
        }
        while (0);

        return error;
    }

    #endif // DOXYGEN_SHOULD_SKIP_THIS


    /******************************************************************//**
     * \name Single-channel samples
     *********************************************************************/
    //@{


    /**
     * \brief Computes a device-wide histogram.  Uses fast block-sorting to compute the histogram. Delivers consistent throughput regardless of sample diversity, but occupancy may be limited by histogram bin count.
     *
     * However, because histograms are privatized in shared memory, a large
     * number of bins (e.g., thousands) may adversely affect occupancy and
     * performance (or even the ability to launch).
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the computation of a 256-bin histogram of
     * single-channel <tt>unsigned char</tt> samples.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input samples and 256-bin output histogram
     * unsigned char *d_samples;
     * unsigned int *d_histogram;
     * int num_items = ...
     * ...
     *
     * // Wrap d_samples device pointer in a random-access texture iterator
     * cub::TexIteratorRA<unsigned int> d_samples_tex_itr;
     * d_samples_tex_itr.BindTexture(d_samples, num_items * sizeof(unsigned char));
     *
     * // Determine temporary device storage requirements for histogram computation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceHistogram::SingleChannelSorting<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histogram, num_items);
     *
     * // Allocate temporary storage for histogram computation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Compute histogram
     * cub::DeviceHistogram::SingleChannelSorting<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histogram, num_items);
     *
     * // Unbind texture iterator
     * d_samples_tex_itr.UnbindTexture();
     *
     * \endcode
     *
     * \tparam BINS                 Number of histogram bins per channel
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)  Must have a value type that can be cast as an integer in the range [0..BINS-1]
     * \tparam HistoCounter         <b>[inferred]</b> Integral type for counting sample occurrences per histogram bin
     */
    template <
        int                 BINS,
        typename            InputIteratorRA,
        typename            HistoCounter>
    __host__ __device__ __forceinline__
    static cudaError_t SingleChannelSorting(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_samples,                          ///< [in] Input samples
        HistoCounter*       d_histogram,                        ///< [out] Array of BINS counters of integral type \p HistoCounter.
        int                 num_samples,                        ///< [in] Number of samples to process
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
        return Dispatch<GRID_HISTO_SORT, BINS, 1, 1>(
            d_temp_storage, temp_storage_bytes, d_samples, &d_histogram, num_samples, stream, stream_synchronous);
    }


    /**
     * \brief Computes a device-wide histogram.  Uses shared-memory atomic read-modify-write operations to compute the histogram.  Input samples having lower diversity can cause performance to be degraded, and occupancy may be limited by histogram bin count.
     *
     * However, because histograms are privatized in shared memory, a large
     * number of bins (e.g., thousands) may adversely affect occupancy and
     * performance (or even the ability to launch).
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the computation of a 256-bin histogram of
     * single-channel <tt>unsigned char</tt> samples.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input samples and 256-bin output histogram
     * unsigned char *d_samples;
     * unsigned int *d_histogram;
     * int num_items = ...
     * ...
     *
     * // Wrap d_samples device pointer in a random-access texture iterator
     * cub::TexIteratorRA<unsigned int> d_samples_tex_itr;
     * d_samples_tex_itr.BindTexture(d_samples, num_items * sizeof(unsigned char));
     *
     * // Determine temporary device storage requirements for histogram computation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceHistogram::SingleChannelSorting<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histogram, num_items);
     *
     * // Allocate temporary storage for histogram computation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Compute histogram
     * cub::DeviceHistogram::SingleChannelSharedAtomic<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histogram, num_items);
     *
     * // Unbind texture iterator
     * d_samples_tex_itr.UnbindTexture();
     *
     * \endcode
     *
     * \tparam BINS                 Number of histogram bins per channel
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)  Must have a value type that can be cast as an integer in the range [0..BINS-1]
     * \tparam HistoCounter         <b>[inferred]</b> Integral type for counting sample occurrences per histogram bin
     */
    template <
        int                 BINS,
        typename            InputIteratorRA,
        typename            HistoCounter>
    __host__ __device__ __forceinline__
    static cudaError_t SingleChannelSharedAtomic(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_samples,                          ///< [in] Input samples
        HistoCounter*       d_histogram,                        ///< [out] Array of BINS counters of integral type \p HistoCounter.
        int                 num_samples,                        ///< [in] Number of samples to process
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch<GRID_HISTO_SHARED_ATOMIC, BINS, 1, 1>(
            d_temp_storage, temp_storage_bytes, d_samples, &d_histogram, num_samples, stream, stream_synchronous);
    }


    /**
     * \brief Computes a device-wide histogram.  Uses global-memory atomic read-modify-write operations to compute the histogram.  Input samples having lower diversity can cause performance to be degraded.
     *
     * Performance is not significantly impacted when computing histograms having large
     * numbers of bins (e.g., thousands).
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the computation of a 256-bin histogram of
     * single-channel <tt>unsigned char</tt> samples.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input samples and 256-bin output histogram
     * unsigned char *d_samples;
     * unsigned int *d_histogram;
     * int num_items = ...
     * ...
     *
     * // Wrap d_samples device pointer in a random-access texture iterator
     * cub::TexIteratorRA<unsigned int> d_samples_tex_itr;
     * d_samples_tex_itr.BindTexture(d_samples, num_items * sizeof(unsigned char));
     *
     * // Determine temporary device storage requirements for histogram computation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceHistogram::SingleChannelSorting<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histogram, num_items);
     *
     * // Allocate temporary storage for histogram computation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Compute histogram
     * cub::DeviceHistogram::SingleChannelGlobalAtomic<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histogram, num_items);
     *
     * // Unbind texture iterator
     * d_samples_tex_itr.UnbindTexture();
     *
     * \endcode
     *
     * \tparam BINS                 Number of histogram bins per channel
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)  Must have a value type that can be cast as an integer in the range [0..BINS-1]
     * \tparam HistoCounter         <b>[inferred]</b> Integral type for counting sample occurrences per histogram bin
     */
    template <
        int                 BINS,
        typename            InputIteratorRA,
        typename            HistoCounter>
    __host__ __device__ __forceinline__
    static cudaError_t SingleChannelGlobalAtomic(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_samples,                          ///< [in] Input samples
        HistoCounter*       d_histogram,                        ///< [out] Array of BINS counters of integral type \p HistoCounter.
        int                 num_samples,                        ///< [in] Number of samples to process
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch<GRID_HISTO_GLOBAL_ATOMIC, BINS, 1, 1>(
            d_temp_storage, temp_storage_bytes, d_samples, &d_histogram, num_samples, stream, stream_synchronous);
    }


    //@}  end member group
    /******************************************************************//**
     * \name Interleaved multi-channel samples
     *********************************************************************/
    //@{


    /**
     * \brief Computes a device-wide histogram from multi-channel data.  Uses fast block-sorting to compute the histogram.  Delivers consistent throughput regardless of sample diversity, but occupancy may be limited by histogram bin count.
     *
     * However, because histograms are privatized in shared memory, a large
     * number of bins (e.g., thousands) may adversely affect occupancy and
     * performance (or even the ability to launch).
     *
     * The total number of samples across all channels (\p num_samples) must be a whole multiple of \p CHANNELS.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the computation of three 256-bin histograms from
     * interleaved quad-channel <tt>unsigned char</tt> samples (e.g., RGB histograms from RGBA samples).
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input samples and
     * // three 256-bin output histograms
     * unsigned char *d_samples;
     * unsigned int *d_histograms[3];
     * int num_items = ...
     * ...
     *
     * // Wrap d_samples device pointer in a random-access texture iterator
     * cub::TexIteratorRA<unsigned int> d_samples_tex_itr;
     * d_samples_tex_itr.BindTexture(d_samples, num_items * sizeof(unsigned char));
     *
     * // Determine temporary device storage requirements for histogram computation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceHistogram::MultiChannelSorting<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histograms, num_items);
     *
     * // Allocate temporary storage for histogram computation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Compute histograms
     * cub::DeviceHistogram::MultiChannelSorting<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histograms, num_items);
     *
     * // Unbind texture iterator
     * d_samples_tex_itr.UnbindTexture();
     *
     * \endcode
     *
     * \tparam BINS                 Number of histogram bins per channel
     * \tparam CHANNELS             Number of channels interleaved in the input data (may be greater than the number of channels being actively histogrammed)
     * \tparam ACTIVE_CHANNELS      <b>[inferred]</b> Number of channels actively being histogrammed
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)  Must have a value type that can be cast as an integer in the range [0..BINS-1]
     * \tparam HistoCounter         <b>[inferred]</b> Integral type for counting sample occurrences per histogram bin
     */
    template <
        int                 BINS,
        int                 CHANNELS,
        int                 ACTIVE_CHANNELS,
        typename            InputIteratorRA,
        typename            HistoCounter>
    __host__ __device__ __forceinline__
    static cudaError_t MultiChannelSorting(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_samples,                          ///< [in] Input samples. The samples from different channels are assumed to be interleaved (e.g., an array of 32b pixels where each pixel consists of four RGBA 8b samples).
        HistoCounter        *d_histograms[ACTIVE_CHANNELS],     ///< [out] Array of channel histogram counter arrays, each having BINS counters of integral type \p HistoCounter.
        int                 num_samples,                        ///< [in] Total number of samples to process in all channels, including non-active channels
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch<GRID_HISTO_SORT, BINS, CHANNELS, ACTIVE_CHANNELS>(
            d_temp_storage, temp_storage_bytes, d_samples, d_histograms, num_samples, stream, stream_synchronous);
    }


    /**
     * \brief Computes a device-wide histogram from multi-channel data.  Uses shared-memory atomic read-modify-write operations to compute the histogram.  Input samples having lower diversity can cause performance to be degraded, and occupancy may be limited by histogram bin count.
     *
     * However, because histograms are privatized in shared memory, a large
     * number of bins (e.g., thousands) may adversely affect occupancy and
     * performance (or even the ability to launch).
     *
     * The total number of samples across all channels (\p num_samples) must be a whole multiple of \p CHANNELS.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the computation of three 256-bin histograms from
     * interleaved quad-channel <tt>unsigned char</tt> samples (e.g., RGB histograms from RGBA samples).
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input samples and
     * // three 256-bin output histograms
     * unsigned char *d_samples;
     * unsigned int *d_histograms[3];
     * int num_items = ...
     * ...
     *
     * // Wrap d_samples device pointer in a random-access texture iterator
     * cub::TexIteratorRA<unsigned int> d_samples_tex_itr;
     * d_samples_tex_itr.BindTexture(d_samples, num_items * sizeof(unsigned char));
     *
     * // Determine temporary device storage requirements for histogram computation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceHistogram::MultiChannelSharedAtomic<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histograms, num_items);
     *
     * // Allocate temporary storage for histogram computation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Compute histograms
     * cub::DeviceHistogram::MultiChannelSharedAtomic<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histograms, num_items);
     *
     * // Unbind texture iterator
     * d_samples_tex_itr.UnbindTexture();
     *
     * \endcode
     *
     * \tparam BINS                 Number of histogram bins per channel
     * \tparam CHANNELS             Number of channels interleaved in the input data (may be greater than the number of channels being actively histogrammed)
     * \tparam ACTIVE_CHANNELS      <b>[inferred]</b> Number of channels actively being histogrammed
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)  Must have a value type that can be cast as an integer in the range [0..BINS-1]
     * \tparam HistoCounter         <b>[inferred]</b> Integral type for counting sample occurrences per histogram bin
     */
    template <
        int                 BINS,
        int                 CHANNELS,
        int                 ACTIVE_CHANNELS,
        typename            InputIteratorRA,
        typename            HistoCounter>
    __host__ __device__ __forceinline__
    static cudaError_t MultiChannelSharedAtomic(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_samples,                          ///< [in] Input samples. The samples from different channels are assumed to be interleaved (e.g., an array of 32b pixels where each pixel consists of four RGBA 8b samples).
        HistoCounter        *d_histograms[ACTIVE_CHANNELS],     ///< [out] Array of channel histogram counter arrays, each having BINS counters of integral type \p HistoCounter.
        int                 num_samples,                        ///< [in] Total number of samples to process in all channels, including non-active channels
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch<GRID_HISTO_SHARED_ATOMIC, BINS, CHANNELS, ACTIVE_CHANNELS>(
            d_temp_storage, temp_storage_bytes, d_samples, d_histograms, num_samples, stream, stream_synchronous);
    }


    /**
     * \brief Computes a device-wide histogram from multi-channel data.  Uses global-memory atomic read-modify-write operations to compute the histogram.  Input samples having lower diversity can cause performance to be degraded.
     *
     * Performance is not significantly impacted when computing histograms having large
     * numbers of bins (e.g., thousands).
     *
     * The total number of samples across all channels (\p num_samples) must be a whole multiple of \p CHANNELS.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * Performance is often improved when referencing input samples through a texture-caching iterator, e.g., cub::TexIteratorRA or cub::TexTransformIteratorRA.
     *
     * \par
     * The code snippet below illustrates the computation of three 256-bin histograms from
     * interleaved quad-channel <tt>unsigned char</tt> samples (e.g., RGB histograms from RGBA samples).
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input samples and
     * // three 256-bin output histograms
     * unsigned char *d_samples;
     * unsigned int *d_histograms[3];
     * int num_items = ...
     * ...
     *
     * // Wrap d_samples device pointer in a random-access texture iterator
     * cub::TexIteratorRA<unsigned int> d_samples_tex_itr;
     * d_samples_tex_itr.BindTexture(d_samples, num_items * sizeof(unsigned char));
     *
     * // Determine temporary device storage requirements for histogram computation
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceHistogram::MultiChannelGlobalAtomic<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histograms, num_items);
     *
     * // Allocate temporary storage for histogram computation
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Compute histograms
     * cub::DeviceHistogram::MultiChannelGlobalAtomic<256>(d_temp_storage, temp_storage_bytes, d_samples_tex_itr, d_histograms, num_items);
     *
     * // Unbind texture iterator
     * d_samples_tex_itr.UnbindTexture();
     *
     * \endcode
     *
     * \tparam BINS                 Number of histogram bins per channel
     * \tparam CHANNELS             Number of channels interleaved in the input data (may be greater than the number of channels being actively histogrammed)
     * \tparam ACTIVE_CHANNELS      <b>[inferred]</b> Number of channels actively being histogrammed
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)  Must have a value type that can be cast as an integer in the range [0..BINS-1]
     * \tparam HistoCounter         <b>[inferred]</b> Integral type for counting sample occurrences per histogram bin
     */
    template <
        int                 BINS,
        int                 CHANNELS,
        int                 ACTIVE_CHANNELS,
        typename            InputIteratorRA,
        typename            HistoCounter>
    __host__ __device__ __forceinline__
    static cudaError_t MultiChannelGlobalAtomic(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_samples,                          ///< [in] Input samples. The samples from different channels are assumed to be interleaved (e.g., an array of 32b pixels where each pixel consists of four RGBA 8b samples).
        HistoCounter        *d_histograms[ACTIVE_CHANNELS],     ///< [out] Array of channel histogram counter arrays, each having BINS counters of integral type \p HistoCounter.
        int                 num_samples,                        ///< [in] Total number of samples to process in all channels, including non-active channels
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch<GRID_HISTO_GLOBAL_ATOMIC, BINS, CHANNELS, ACTIVE_CHANNELS>(
            d_temp_storage, temp_storage_bytes, d_samples, d_histograms, num_samples, stream, stream_synchronous);
    }

    //@}  end member group

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)


