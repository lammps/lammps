
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
 * cub::DeviceReorder provides device-wide operations for partitioning and filtering lists of items residing within global memory.
 */

#pragma once

#include <stdio.h>
#include <iterator>

#include "device_scan.cuh"
#include "block/block_partition_tiles.cuh"
#include "../grid/grid_queue.cuh"
#include "../util_debug.cuh"
#include "../util_device.cuh"
#include "../util_vector.cuh"
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
 * Partition kernel entry point (multi-block)
 */
template <
    typename    BlockPartitionTilesPolicy,  ///< Tuning policy for cub::BlockPartitionTiles abstraction
    typename    InputIteratorRA,            ///< Random-access iterator type for input (may be a simple pointer type)
    typename    OutputIteratorRA,           ///< Random-access iterator type for output (may be a simple pointer type)
    typename    LengthOutputIterator,       ///< Output iterator type for recording the length of the first partition (may be a simple pointer type)
    typename    PredicateOp,                ///< Unary predicate operator indicating membership in the first partition type having member <tt>bool operator()(const T &val)</tt>
    typename    SizeT>                      ///< Integer type used for global array indexing
__launch_bounds__ (int(BlockPartitionTilesPolicy::BLOCK_THREADS))
__global__ void PartitionKernel(
    InputIteratorRA                                                                         d_in,               ///< Input data
    OutputIteratorRA                                                                        d_out,              ///< Output data
    LengthOutputIterator                                                                    d_partition_length, ///< Number of items in the first partition
    ScanTileDescriptor<PartitionScanTuple<SizeT, BlockPartitionTilesPolicy::PARTITOINS> >   *d_tile_status,     ///< Global list of tile status
    PredicateOp                                                                             pred_op,            ///< Unary predicate operator indicating membership in the first partition
    SizeT                                                                                   num_items,          ///< Total number of input items for the entire problem
    int                                                                                     num_tiles,          ///< Totla number of intut tiles for the entire problem
    GridQueue<int>                                                                          queue)              ///< Descriptor for performing dynamic mapping of tile data to thread blocks
{
    enum
    {
        TILE_STATUS_PADDING = PtxArchProps::WARP_THREADS,
    };

    typedef PartitionScanTuple<SizeT, BlockPartitionTilesPolicy::PARTITOINS> PartitionScanTuple;

    // Thread block type for scanning input tiles
    typedef BlockPartitionTiles<
        BlockPartitionTilesPolicy,
        InputIteratorRA,
        OutputIteratorRA,
        PredicateOp,
        SizeT> BlockPartitionTilesT;

    // Shared memory for BlockPartitionTiles
    __shared__ typename BlockPartitionTilesT::TempStorage temp_storage;

    // Process tiles
    PartitionScanTuple  partition_ends;     // Ending offsets for partitions (one-after)
    bool                is_last_tile;       // Whether or not this block handled the last tile (i.e., partition_ends is valid for the entire input)
    BlockPartitionTilesT(temp_storage, d_in, d_out, d_tile_status + TILE_STATUS_PADDING, pred_op, num_items).ConsumeTiles(
        queue,
        num_tiles,
        partition_ends,
        is_last_tile);

    // Record the length of the first partition
    if (is_last_tile && (threadIdx.x == 0))
    {
        *d_partition_length = partition_ends.x;
    }
}


#endif // DOXYGEN_SHOULD_SKIP_THIS



/******************************************************************************
 * DeviceReorder
 *****************************************************************************/

/**
 * \addtogroup DeviceModule
 * @{
 */

/**
 * \brief DeviceReorder provides device-wide operations for partitioning and filtering lists of items residing within global memory
 */
struct DeviceReorder
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    /// Generic structure for encapsulating dispatch properties.  Mirrors the constants within BlockPartitionTilesPolicy.
    struct KernelDispachParams
    {
        int                     block_threads;
        int                     items_per_thread;
        BlockScanAlgorithm      scan_algorithm;
        int                     tile_size;

        template <typename BlockPartitionTilesPolicy>
        __host__ __device__ __forceinline__
        void Init()
        {
            block_threads               = BlockPartitionTilesPolicy::BLOCK_THREADS;
            items_per_thread            = BlockPartitionTilesPolicy::ITEMS_PER_THREAD;
            scan_algorithm              = BlockPartitionTilesPolicy::SCAN_ALGORITHM;
            tile_size                   = block_threads * items_per_thread;
        }
    };


    /******************************************************************************
     * Tuning policies
     ******************************************************************************/


    /// Specializations of tuned policy types for different PTX architectures
    template <
        int         PARTITIONS,
        typename    T,
        typename    SizeT,
        int         ARCH>
    struct TunedPolicies;

    /// SM35 tune
    template <int PARTITIONS, typename T, typename SizeT>
    struct TunedPolicies<PARTITIONS, T, SizeT, 350>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 16,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };

        typedef BlockPartitionTilesPolicy<PARTITIONS, 128, ITEMS_PER_THREAD, LOAD_LDG, BLOCK_SCAN_RAKING_MEMOIZE> PartitionPolicy;
    };

    /// SM30 tune
    template <int PARTITIONS, typename T, typename SizeT>
    struct TunedPolicies<PARTITIONS, T, SizeT, 300>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 9,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };

        typedef BlockPartitionTilesPolicy<PARTITIONS, 256, ITEMS_PER_THREAD, LOAD_DEFAULT, BLOCK_SCAN_RAKING_MEMOIZE> PartitionPolicy;
    };

    /// SM20 tune
    template <int PARTITIONS, typename T, typename SizeT>
    struct TunedPolicies<PARTITIONS, T, SizeT, 200>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 15,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };

        typedef BlockPartitionTilesPolicy<PARTITIONS, 128, ITEMS_PER_THREAD, LOAD_DEFAULT, BLOCK_SCAN_RAKING_MEMOIZE> PartitionPolicy;
    };

    /// SM10 tune
    template <int PARTITIONS, typename T, typename SizeT>
    struct TunedPolicies<PARTITIONS, T, SizeT, 100>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 7,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };
        typedef BlockPartitionTilesPolicy<PARTITIONS, 128, ITEMS_PER_THREAD, LOAD_DEFAULT, BLOCK_SCAN_RAKING> PartitionPolicy;
    };


    /// Tuning policy for the PTX architecture that DevicePartition operations will get dispatched to
    template <int PARTITIONS, typename T, typename SizeT>
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
        typedef TunedPolicies<PARTITIONS, T, SizeT, PTX_TUNE_ARCH> PtxTunedPolicies;

        // PartitionPolicy that opaquely derives from the specialization corresponding to the current PTX compiler pass
        struct PartitionPolicy : PtxTunedPolicies::PartitionPolicy {};

        /**
         * Initialize dispatch params with the policies corresponding to the PTX assembly we will use
         */
        static void InitDispatchParams(int ptx_version, KernelDispachParams &scan_dispatch_params)
        {
            if (ptx_version >= 350)
            {
                typedef TunedPolicies<PARTITIONS, T, SizeT, 350> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::PartitionPolicy>();
            }
            else if (ptx_version >= 300)
            {
                typedef TunedPolicies<PARTITIONS, T, SizeT, 300> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::PartitionPolicy>();
            }
            else if (ptx_version >= 200)
            {
                typedef TunedPolicies<PARTITIONS, T, SizeT, 200> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::PartitionPolicy>();
            }
            else
            {
                typedef TunedPolicies<PARTITIONS, T, SizeT, 100> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::PartitionPolicy>();
            }
        }
    };


    /******************************************************************************
     * Utility methods
     ******************************************************************************/

    /**
     * Internal dispatch routine
     */
    template <
        typename                    ScanInitKernelPtr,              ///< Function type of cub::ScanInitKernel
        typename                    PartitionKernelPtr,             ///< Function type of cub::PartitionKernel
        typename                    InputIteratorRA,                ///< Random-access iterator type for input (may be a simple pointer type)
        typename                    OutputIteratorRA,               ///< Random-access iterator type for output (may be a simple pointer type)
        typename                    LengthOutputIterator,           ///< Output iterator type for recording the length of the first partition (may be a simple pointer type)
        typename                    PredicateOp,                    ///< Unary predicate operator indicating membership in the first partition type having member <tt>bool operator()(const T &val)</tt>
        typename                    SizeT>                          ///< Integer type used for global array indexing
    __host__ __device__ __forceinline__
    static cudaError_t Dispatch(
        int                         ptx_version,                    ///< [in] PTX version
        void                        *d_temp_storage,                ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t                      &temp_storage_bytes,            ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        ScanInitKernelPtr           init_kernel,                    ///< [in] Kernel function pointer to parameterization of cub::PartitionInitKernel
        PartitionKernelPtr          partition_kernel,               ///< [in] Kernel function pointer to parameterization of cub::PartitionKernel
        KernelDispachParams         &scan_dispatch_params,          ///< [in] Dispatch parameters that match the policy that \p partition_kernel was compiled for
        InputIteratorRA             d_in,                           ///< [in] Iterator pointing to scan input
        OutputIteratorRA            d_out,                          ///< [in] Iterator pointing to scan output
        LengthOutputIterator        d_partition_length,                 ///< [out] Output iterator referencing the location where the pivot offset (i.e., the length of the first partition) is to be recorded
        PredicateOp                 pred_op,                        ///< [in] Unary predicate operator indicating membership in the first partition
        SizeT                       num_items,                      ///< [in] Total number of items to partition
        cudaStream_t                stream              = 0,        ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                        stream_synchronous  = false)    ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {

#ifndef CUB_RUNTIME_ENABLED

        // Kernel launch not supported from this device
        return CubDebug(cudaErrorNotSupported);

#else

        enum
        {
            TILE_STATUS_PADDING = 32,
        };

        // Data type
        typedef typename std::iterator_traits<InputIteratorRA>::value_type T;

        // Scan tuple type and tile status descriptor type
        typedef typename VectorHelper<SizeT, 2>::Type ScanTuple;
        typedef ScanTileDescriptor<ScanTuple> ScanTileDescriptorT;

        cudaError error = cudaSuccess;
        do
        {
            // Number of input tiles
            int num_tiles = (num_items + scan_dispatch_params.tile_size - 1) / scan_dispatch_params.tile_size;

            // Temporary storage allocation requirements
            void* allocations[2];
            size_t allocation_sizes[2] =
            {
                (num_tiles + TILE_STATUS_PADDING) * sizeof(ScanTileDescriptorT),      // bytes needed for tile status descriptors
                GridQueue<int>::AllocationSize()                                            // bytes needed for grid queue descriptor
            };

            // Alias temporaries (or set the necessary size of the storage allocation)
            if (CubDebug(error = AliasTemporaries(d_temp_storage, temp_storage_bytes, allocations, allocation_sizes))) break;

            // Return if the caller is simply requesting the size of the storage allocation
            if (d_temp_storage == NULL)
                return cudaSuccess;

            // Global list of tile status
            ScanTileDescriptorT *d_tile_status = (ScanTileDescriptorT*) allocations[0];

            // Grid queue descriptor
            GridQueue<int> queue(allocations[1]);

            // Log init_kernel configuration
            int init_kernel_threads = 128;
            int init_grid_size = (num_tiles + init_kernel_threads - 1) / init_kernel_threads;
            if (stream_synchronous) CubLog("Invoking init_kernel<<<%d, %d, 0, %lld>>>()\n", init_grid_size, init_kernel_threads, (long long) stream);

            // Invoke init_kernel to initialize tile descriptors and queue descriptors
            init_kernel<<<init_grid_size, init_kernel_threads, 0, stream>>>(
                queue,
                d_tile_status,
                num_tiles);

            // Sync the stream if specified
            if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;

            // Get grid size for multi-block kernel
            int scan_grid_size;
            int multi_sm_occupancy = -1;
            if (ptx_version < 200)
            {
                // We don't have atomics (or don't have fast ones), so just assign one
                // block per tile (limited to 65K tiles)
                scan_grid_size = num_tiles;
            }
            else
            {
                // We have atomics and can thus reuse blocks across multiple tiles using a queue descriptor.
                // Get GPU id
                int device_ordinal;
                if (CubDebug(error = cudaGetDevice(&device_ordinal))) break;

                // Get SM count
                int sm_count;
                if (CubDebug(error = cudaDeviceGetAttribute (&sm_count, cudaDevAttrMultiProcessorCount, device_ordinal))) break;

                // Get a rough estimate of partition_kernel SM occupancy based upon the maximum SM occupancy of the targeted PTX architecture
                multi_sm_occupancy = CUB_MIN(
                    ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADBLOCKS,
                    ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADS / scan_dispatch_params.block_threads);

#ifndef __CUDA_ARCH__
                // We're on the host, so come up with a
                Device device_props;
                if (CubDebug(error = device_props.Init(device_ordinal))) break;

                if (CubDebug(error = device_props.MaxSmOccupancy(
                    multi_sm_occupancy,
                    partition_kernel,
                    scan_dispatch_params.block_threads))) break;
#endif
                // Get device occupancy for partition_kernel
                int scan_occupancy = multi_sm_occupancy * sm_count;

                // Get grid size for partition_kernel
                scan_grid_size = (num_tiles < scan_occupancy) ?
                    num_tiles :                 // Not enough to fill the device with threadblocks
                    scan_occupancy;      // Fill the device with threadblocks
            }

            // Log partition_kernel configuration
            if (stream_synchronous) CubLog("Invoking partition_kernel<<<%d, %d, 0, %lld>>>(), %d items per thread, %d SM occupancy\n",
                scan_grid_size, scan_dispatch_params.block_threads, (long long) stream, scan_dispatch_params.items_per_thread, multi_sm_occupancy);

            // Invoke partition_kernel
            partition_kernel<<<scan_grid_size, scan_dispatch_params.block_threads, 0, stream>>>(
                d_in,
                d_out,
                d_partition_length,
                d_tile_status,
                pred_op,
                num_items,
                num_tiles,
                queue);

            // Sync the stream if specified
            if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;
        }
        while (0);

        return error;

#endif  // CUB_RUNTIME_ENABLED
    }



    /**
     * Internal partition dispatch routine for using default tuning policies
     */
    template <
        typename                    PARTITIONS,                     ///< Number of partitions we are keeping
        typename                    InputIteratorRA,                ///< Random-access iterator type for input (may be a simple pointer type)
        typename                    OutputIteratorRA,               ///< Random-access iterator type for output (may be a simple pointer type)
        typename                    LengthOutputIterator,           ///< Output iterator type for recording the length of the first partition (may be a simple pointer type)
        typename                    PredicateOp,                    ///< Unary predicate operator indicating membership in the first partition type having member <tt>bool operator()(const T &val)</tt>
        typename                    SizeT>                          ///< Integer type used for global array indexing
    __host__ __device__ __forceinline__
    static cudaError_t Dispatch(
        void                        *d_temp_storage,                ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t                      &temp_storage_bytes,            ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA             d_in,                           ///< [in] Iterator pointing to input items
        OutputIteratorRA            d_out,                          ///< [in] Iterator pointing to output items
        LengthOutputIterator        d_partition_length,             ///< [out] Output iterator referencing the location where the pivot offset (i.e., the length of the first partition) is to be recorded
        PredicateOp                 pred_op,                        ///< [in] Unary predicate operator indicating membership in the first partition
        SizeT                       num_items,                      ///< [in] Total number of items to partition
        cudaStream_t                stream              = 0,        ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                        stream_synchronous  = false)    ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
        // Data type
        typedef typename std::iterator_traits<InputIteratorRA>::value_type T;

        // Tuning polices
        typedef PtxDefaultPolicies<PARTITIONS, T, SizeT>        PtxDefaultPolicies;     // Wrapper of default kernel policies
        typedef typename PtxDefaultPolicies::PartitionPolicy    PartitionPolicy;        // Partition kernel policy

        cudaError error = cudaSuccess;
        do
        {
            // Declare dispatch parameters
            KernelDispachParams scan_dispatch_params;

            int ptx_version;
#ifdef __CUDA_ARCH__
            // We're on the device, so initialize the dispatch parameters with the PtxDefaultPolicies directly
            scan_dispatch_params.Init<PartitionPolicy>();
            ptx_version = CUB_PTX_ARCH;
#else
            // We're on the host, so lookup and initialize the dispatch parameters with the policies that match the device's PTX version
            if (CubDebug(error = PtxVersion(ptx_version))) break;
            PtxDefaultPolicies::InitDispatchParams(ptx_version, scan_dispatch_params);
#endif

            Dispatch(
                ptx_version,
                d_temp_storage,
                temp_storage_bytes,
                ScanInitKernel<T, SizeT>,
                PartitionKernel<PartitionPolicy, InputIteratorRA, OutputIteratorRA, LengthOutputIterator, PredicateOp, SizeT>,
                scan_dispatch_params,
                d_in,
                d_out,
                d_partition_length,
                pred_op,
                num_items,
                stream,
                stream_synchronous);

            if (CubDebug(error)) break;
        }
        while (0);

        return error;
    }

    #endif // DOXYGEN_SHOULD_SKIP_THIS


    /**
     * \brief Splits a list of input items into two partitions within the given output list using the specified predicate.  The relative ordering of inputs is not necessarily preserved.
     *
     * An item \p val is placed in the first partition if <tt>pred_op(val) == true</tt>, otherwise
     * it is placed in the second partition.  The offset of the partitioning pivot (equivalent to
     * the total length of the first partition as well as the starting offset of the second), is
     * recorded to \p d_partition_length.
     *
     * The length of the output referenced by \p d_out is assumed to be the same as that of \p d_in.
     *
     * \devicestorage
     *
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)
     * \tparam OutputIteratorRA     <b>[inferred]</b> Random-access iterator type for output (may be a simple pointer type)
     * \tparam LengthOutputIterator <b>[inferred]</b> Random-access iterator type for output (may be a simple pointer type)
     * \tparam PredicateOp          <b>[inferred]</b> Unary predicate operator indicating membership in the first partition type having member <tt>bool operator()(const T &val)</tt>
     */
    template <
        typename                InputIteratorRA,
        typename                OutputIteratorRA,
        typename                LengthOutputIterator,
        typename                PredicateOp>
    __host__ __device__ __forceinline__
    static cudaError_t Partition(
        void                    *d_temp_storage,                ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t                  &temp_storage_bytes,            ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA         d_in,                           ///< [in] Iterator pointing to input items
        OutputIteratorRA        d_out,                          ///< [in] Iterator pointing to output items
        LengthOutputIterator    d_pivot_offset,                 ///< [out] Output iterator referencing the location where the pivot offset is to be recorded
        PredicateOp             pred_op,                        ///< [in] Unary predicate operator indicating membership in the first partition
        int                     num_items,                      ///< [in] Total number of items to partition
        cudaStream_t            stream              = 0,        ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                    stream_synchronous  = false)    ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        typedef typename std::iterator_traits<InputIteratorRA>::value_type T;
        return Dispatch(d_temp_storage, temp_storage_bytes, d_in, d_out, Sum(), T(), num_items, stream, stream_synchronous);
    }


};


/** @} */       // DeviceModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)


