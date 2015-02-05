
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
 * cub::DeviceScan provides operations for computing a device-wide, parallel prefix scan across data items residing within global memory.
 */

#pragma once

#include <stdio.h>
#include <iterator>

#include "block/block_scan_tiles.cuh"
#include "../thread/thread_operators.cuh"
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
 * Initialization kernel for tile status initialization (multi-block)
 */
template <
    typename T,                                     ///< Scan value type
    typename SizeT>                                 ///< Integer type used for global array indexing
__global__ void ScanInitKernel(
    GridQueue<SizeT>            grid_queue,         ///< [in] Descriptor for performing dynamic mapping of input tiles to thread blocks
    ScanTileDescriptor<T>       *d_tile_status,     ///< [out] Tile status words
    int                         num_tiles)          ///< [in] Number of tiles
{
    typedef ScanTileDescriptor<T> ScanTileDescriptorT;

    enum
    {
        TILE_STATUS_PADDING = PtxArchProps::WARP_THREADS,
    };

    // Reset queue descriptor
    if ((blockIdx.x == 0) && (threadIdx.x == 0)) grid_queue.ResetDrain(num_tiles);

    // Initialize tile status
    int tile_offset = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (tile_offset < num_tiles)
    {
        // Not-yet-set
        d_tile_status[TILE_STATUS_PADDING + tile_offset].status = SCAN_TILE_INVALID;
    }

    if ((blockIdx.x == 0) && (threadIdx.x < TILE_STATUS_PADDING))
    {
        // Padding
        d_tile_status[threadIdx.x].status = SCAN_TILE_OOB;
    }
}


/**
 * Scan kernel entry point (multi-block)
 */
template <
    typename    BlockScanTilesPolicy,           ///< Tuning policy for cub::BlockScanTiles abstraction
    typename    InputIteratorRA,                ///< Random-access iterator type for input (may be a simple pointer type)
    typename    OutputIteratorRA,               ///< Random-access iterator type for output (may be a simple pointer type)
    typename    T,                              ///< The scan data type
    typename    ScanOp,                         ///< Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
    typename    Identity,                       ///< Identity value type (cub::NullType for inclusive scans)
    typename    SizeT>                          ///< Integer type used for global array indexing
__launch_bounds__ (int(BlockScanTilesPolicy::BLOCK_THREADS))
__global__ void ScanKernel(
    InputIteratorRA             d_in,           ///< Input data
    OutputIteratorRA            d_out,          ///< Output data
    ScanTileDescriptor<T>       *d_tile_status, ///< Global list of tile status
    ScanOp                      scan_op,        ///< Binary scan operator
    Identity                    identity,       ///< Identity element
    SizeT                       num_items,      ///< Total number of scan items for the entire problem
    GridQueue<int>              queue)          ///< Descriptor for performing dynamic mapping of tile data to thread blocks
{
    enum
    {
        TILE_STATUS_PADDING = PtxArchProps::WARP_THREADS,
    };

    // Thread block type for scanning input tiles
    typedef BlockScanTiles<
        BlockScanTilesPolicy,
        InputIteratorRA,
        OutputIteratorRA,
        ScanOp,
        Identity,
        SizeT> BlockScanTilesT;

    // Shared memory for BlockScanTiles
    __shared__ typename BlockScanTilesT::TempStorage temp_storage;

    // Process tiles
    BlockScanTilesT(temp_storage, d_in, d_out, scan_op, identity).ConsumeTiles(
        num_items,
        queue,
        d_tile_status + TILE_STATUS_PADDING);
}


#endif // DOXYGEN_SHOULD_SKIP_THIS



/******************************************************************************
 * DeviceScan
 *****************************************************************************/

/**
 * \brief DeviceScan provides operations for computing a device-wide, parallel prefix scan across data items residing within global memory. ![](device_scan.png)
 * \ingroup DeviceModule
 *
 * \par Overview
 * Given a list of input elements and a binary reduction operator, a [<em>prefix scan</em>](http://en.wikipedia.org/wiki/Prefix_sum)
 * produces an output list where each element is computed to be the reduction
 * of the elements occurring earlier in the input list.  <em>Prefix sum</em>
 * connotes a prefix scan with the addition operator. The term \em inclusive indicates
 * that the <em>i</em><sup>th</sup> output reduction incorporates the <em>i</em><sup>th</sup> input.
 * The term \em exclusive indicates the <em>i</em><sup>th</sup> input is not incorporated into
 * the <em>i</em><sup>th</sup> output reduction.
 *
 * \par Usage Considerations
 * \cdp_class{DeviceScan}
 *
 * \par Performance
 *
 * \image html scan_perf.png
 *
 */
struct DeviceScan
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    /// Generic structure for encapsulating dispatch properties.  Mirrors the constants within BlockScanTilesPolicy.
    struct KernelDispachParams
    {
        // Policy fields
        int                     block_threads;
        int                     items_per_thread;
        BlockLoadAlgorithm      load_policy;
        BlockStoreAlgorithm     store_policy;
        BlockScanAlgorithm      scan_algorithm;

        // Other misc
        int                     tile_size;

        template <typename BlockScanTilesPolicy>
        __host__ __device__ __forceinline__
        void Init()
        {
            block_threads               = BlockScanTilesPolicy::BLOCK_THREADS;
            items_per_thread            = BlockScanTilesPolicy::ITEMS_PER_THREAD;
            load_policy                 = BlockScanTilesPolicy::LOAD_ALGORITHM;
            store_policy                = BlockScanTilesPolicy::STORE_ALGORITHM;
            scan_algorithm              = BlockScanTilesPolicy::SCAN_ALGORITHM;

            tile_size                   = block_threads * items_per_thread;
        }

        __host__ __device__ __forceinline__
        void Print()
        {
            printf("%d, %d, %d, %d, %d",
                block_threads,
                items_per_thread,
                load_policy,
                store_policy,
                scan_algorithm);
        }

    };


    /******************************************************************************
     * Tuning policies
     ******************************************************************************/


    /// Specializations of tuned policy types for different PTX architectures
    template <
        typename    T,
        typename    SizeT,
        int         ARCH>
    struct TunedPolicies;

    /// SM35 tune
    template <typename T, typename SizeT>
    struct TunedPolicies<T, SizeT, 350>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 16,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };

        // ScanPolicy: GTX Titan: 29.1B items/s (232.4 GB/s) @ 48M 32-bit T
        typedef BlockScanTilesPolicy<128, ITEMS_PER_THREAD,  BLOCK_LOAD_DIRECT, false, LOAD_LDG, BLOCK_STORE_WARP_TRANSPOSE, true, BLOCK_SCAN_RAKING_MEMOIZE> ScanPolicy;
    };

    /// SM30 tune
    template <typename T, typename SizeT>
    struct TunedPolicies<T, SizeT, 300>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 9,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };

        typedef BlockScanTilesPolicy<256, ITEMS_PER_THREAD,  BLOCK_LOAD_WARP_TRANSPOSE, false, LOAD_DEFAULT, BLOCK_STORE_WARP_TRANSPOSE, false, BLOCK_SCAN_RAKING_MEMOIZE> ScanPolicy;
    };

    /// SM20 tune
    template <typename T, typename SizeT>
    struct TunedPolicies<T, SizeT, 200>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 15,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };

        // ScanPolicy: GTX 580: 20.3B items/s (162.3 GB/s) @ 48M 32-bit T
        typedef BlockScanTilesPolicy<128, ITEMS_PER_THREAD, BLOCK_LOAD_WARP_TRANSPOSE, false, LOAD_DEFAULT, BLOCK_STORE_WARP_TRANSPOSE, false, BLOCK_SCAN_RAKING_MEMOIZE> ScanPolicy;
    };

    /// SM10 tune
    template <typename T, typename SizeT>
    struct TunedPolicies<T, SizeT, 100>
    {
        enum {
            NOMINAL_4B_ITEMS_PER_THREAD = 7,
            ITEMS_PER_THREAD            = CUB_MIN(NOMINAL_4B_ITEMS_PER_THREAD, CUB_MAX(1, (NOMINAL_4B_ITEMS_PER_THREAD * 4 / sizeof(T)))),
        };
        typedef BlockScanTilesPolicy<128, ITEMS_PER_THREAD, BLOCK_LOAD_TRANSPOSE, false, LOAD_DEFAULT, BLOCK_STORE_TRANSPOSE, false, BLOCK_SCAN_RAKING> ScanPolicy;
    };


    /// Tuning policy for the PTX architecture that DeviceScan operations will get dispatched to
    template <typename T, typename SizeT>
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
        typedef TunedPolicies<T, SizeT, PTX_TUNE_ARCH> PtxTunedPolicies;

        // ScanPolicy that opaquely derives from the specialization corresponding to the current PTX compiler pass
        struct ScanPolicy : PtxTunedPolicies::ScanPolicy {};

        /**
         * Initialize dispatch params with the policies corresponding to the PTX assembly we will use
         */
        static void InitDispatchParams(int ptx_version, KernelDispachParams &scan_dispatch_params)
        {
            if (ptx_version >= 350)
            {
                typedef TunedPolicies<T, SizeT, 350> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::ScanPolicy>();
            }
            else if (ptx_version >= 300)
            {
                typedef TunedPolicies<T, SizeT, 300> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::ScanPolicy>();
            }
            else if (ptx_version >= 200)
            {
                typedef TunedPolicies<T, SizeT, 200> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::ScanPolicy>();
            }
            else
            {
                typedef TunedPolicies<T, SizeT, 100> TunedPolicies;
                scan_dispatch_params.Init<typename TunedPolicies::ScanPolicy>();
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
        typename                    ScanKernelPtr,                  ///< Function type of cub::ScanKernel
        typename                    InputIteratorRA,                ///< Random-access iterator type for input (may be a simple pointer type)
        typename                    OutputIteratorRA,               ///< Random-access iterator type for output (may be a simple pointer type)
        typename                    ScanOp,                         ///< Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
        typename                    Identity,                       ///< Identity value type (cub::NullType for inclusive scans)
        typename                    SizeT>                          ///< Integer type used for global array indexing
    __host__ __device__ __forceinline__
    static cudaError_t Dispatch(
        int                         ptx_version,                    ///< [in] PTX version
        void                        *d_temp_storage,                ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t                      &temp_storage_bytes,            ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        ScanInitKernelPtr           init_kernel,                    ///< [in] Kernel function pointer to parameterization of cub::ScanInitKernel
        ScanKernelPtr               scan_kernel,                    ///< [in] Kernel function pointer to parameterization of cub::ScanKernel
        KernelDispachParams         &scan_dispatch_params,          ///< [in] Dispatch parameters that match the policy that \p scan_kernel was compiled for
        InputIteratorRA             d_in,                           ///< [in] Iterator pointing to scan input
        OutputIteratorRA            d_out,                          ///< [in] Iterator pointing to scan output
        ScanOp                      scan_op,                        ///< [in] Binary scan operator
        Identity                    identity,                       ///< [in] Identity element
        SizeT                       num_items,                      ///< [in] Total number of items to scan
        cudaStream_t                stream              = 0,        ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                        stream_synchronous  = false)    ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {

#ifndef CUB_RUNTIME_ENABLED

        // Kernel launch not supported from this device
        return CubDebug(cudaErrorNotSupported);

#else

        enum
        {
            TILE_STATUS_PADDING     = 32,
            INIT_KERNEL_THREADS     = 128
        };

        // Data type
        typedef typename std::iterator_traits<InputIteratorRA>::value_type T;

        // Tile status descriptor type
        typedef ScanTileDescriptor<T> ScanTileDescriptorT;

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
                GridQueue<int>::AllocationSize()                                      // bytes needed for grid queue descriptor
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
            int init_grid_size = (num_tiles + INIT_KERNEL_THREADS - 1) / INIT_KERNEL_THREADS;
            if (stream_synchronous) CubLog("Invoking init_kernel<<<%d, %d, 0, %lld>>>()\n", init_grid_size, INIT_KERNEL_THREADS, (long long) stream);

            // Invoke init_kernel to initialize tile descriptors and queue descriptors
            init_kernel<<<init_grid_size, INIT_KERNEL_THREADS, 0, stream>>>(
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

                // Get a rough estimate of scan_kernel SM occupancy based upon the maximum SM occupancy of the targeted PTX architecture
                multi_sm_occupancy = CUB_MIN(
                    ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADBLOCKS,
                    ArchProps<CUB_PTX_ARCH>::MAX_SM_THREADS / scan_dispatch_params.block_threads);

#ifndef __CUDA_ARCH__
                // We're on the host, so come up with a
                Device device_props;
                if (CubDebug(error = device_props.Init(device_ordinal))) break;

                if (CubDebug(error = device_props.MaxSmOccupancy(
                    multi_sm_occupancy,
                    scan_kernel,
                    scan_dispatch_params.block_threads))) break;
#endif
                // Get device occupancy for scan_kernel
                int scan_occupancy = multi_sm_occupancy * sm_count;

                // Get grid size for scan_kernel
                scan_grid_size = (num_tiles < scan_occupancy) ?
                    num_tiles :                 // Not enough to fill the device with threadblocks
                    scan_occupancy;      // Fill the device with threadblocks
            }

            // Log scan_kernel configuration
            if (stream_synchronous) CubLog("Invoking scan_kernel<<<%d, %d, 0, %lld>>>(), %d items per thread, %d SM occupancy\n",
                scan_grid_size, scan_dispatch_params.block_threads, (long long) stream, scan_dispatch_params.items_per_thread, multi_sm_occupancy);

            // Invoke scan_kernel
            scan_kernel<<<scan_grid_size, scan_dispatch_params.block_threads, 0, stream>>>(
                d_in,
                d_out,
                d_tile_status,
                scan_op,
                identity,
                num_items,
                queue);

            // Sync the stream if specified
            if (stream_synchronous && (CubDebug(error = SyncStream(stream)))) break;
        }
        while (0);

        return error;

#endif  // CUB_RUNTIME_ENABLED
    }



    /**
     * Internal scan dispatch routine for using default tuning policies
     */
    template <
        typename                    InputIteratorRA,                ///< Random-access iterator type for input (may be a simple pointer type)
        typename                    OutputIteratorRA,               ///< Random-access iterator type for output (may be a simple pointer type)
        typename                    ScanOp,                         ///< Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
        typename                    Identity,                       ///< Identity value type (cub::NullType for inclusive scans)
        typename                    SizeT>                          ///< Integer type used for global array indexing
    __host__ __device__ __forceinline__
    static cudaError_t Dispatch(
        void                        *d_temp_storage,                ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t                      &temp_storage_bytes,            ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA             d_in,                           ///< [in] Iterator pointing to scan input
        OutputIteratorRA            d_out,                          ///< [in] Iterator pointing to scan output
        ScanOp                      scan_op,                        ///< [in] Binary scan operator
        Identity                    identity,                       ///< [in] Identity element
        SizeT                       num_items,                      ///< [in] Total number of items to scan
        cudaStream_t                stream              = 0,        ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                        stream_synchronous  = false)    ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Default is \p false.
    {
        // Data type
        typedef typename std::iterator_traits<InputIteratorRA>::value_type T;

        // Tuning polices
        typedef PtxDefaultPolicies<T, SizeT>                    PtxDefaultPolicies;     // Wrapper of default kernel policies
        typedef typename PtxDefaultPolicies::ScanPolicy   ScanPolicy;       // Scan kernel policy

        cudaError error = cudaSuccess;
        do
        {
            // Declare dispatch parameters
            KernelDispachParams scan_dispatch_params;

            int ptx_version;
#ifdef __CUDA_ARCH__
            // We're on the device, so initialize the dispatch parameters with the PtxDefaultPolicies directly
            scan_dispatch_params.Init<ScanPolicy>();
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
                ScanKernel<ScanPolicy, InputIteratorRA, OutputIteratorRA, T, ScanOp, Identity, SizeT>,
                scan_dispatch_params,
                d_in,
                d_out,
                scan_op,
                identity,
                num_items,
                stream,
                stream_synchronous);

            if (CubDebug(error)) break;
        }
        while (0);

        return error;
    }

    #endif // DOXYGEN_SHOULD_SKIP_THIS


    /******************************************************************//**
     * \name Exclusive scans
     *********************************************************************/
    //@{

    /**
     * \brief Computes a device-wide exclusive prefix sum.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the exclusive prefix sum of a device vector of \p int items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input and output
     * int *d_scan_input, *d_scan_output;
     * int num_items = ...
     *
     * ...
     *
     * // Determine temporary device storage requirements for exclusive prefix sum
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, num_items);
     *
     * // Allocate temporary storage for exclusive prefix sum
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Run exclusive prefix sum
     * cub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, num_items);
     *
     * \endcode
     *
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)
     * \tparam OutputIteratorRA     <b>[inferred]</b> Random-access iterator type for output (may be a simple pointer type)
     */
    template <
        typename            InputIteratorRA,
        typename            OutputIteratorRA>
    __host__ __device__ __forceinline__
    static cudaError_t ExclusiveSum(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_in,                               ///< [in] Iterator pointing to scan input
        OutputIteratorRA    d_out,                              ///< [in] Iterator pointing to scan output
        int                 num_items,                          ///< [in] Total number of items to scan
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        typedef typename std::iterator_traits<InputIteratorRA>::value_type T;
        return Dispatch(d_temp_storage, temp_storage_bytes, d_in, d_out, Sum(), T(), num_items, stream, stream_synchronous);
    }


    /**
     * \brief Computes a device-wide exclusive prefix scan using the specified binary \p scan_op functor.
     *
     * \par
     * Supports non-commutative scan operators.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the exclusive prefix scan of a device vector of \p int items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input and output
     * int *d_scan_input, *d_scan_output;
     * int num_items = ...
     *
     * ...
     *
     * // Determine temporary device storage requirements for exclusive prefix scan
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceScan::ExclusiveScan(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, cub::Max(), (int) MIN_INT, num_items);
     *
     * // Allocate temporary storage for exclusive prefix scan
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Run exclusive prefix scan (max)
     * cub::DeviceScan::ExclusiveScan(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, cub::Max(), (int) MIN_INT, num_items);
     *
     * \endcode
     *
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)
     * \tparam OutputIteratorRA     <b>[inferred]</b> Random-access iterator type for output (may be a simple pointer type)
     * \tparam ScanOp               <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     * \tparam Identity             <b>[inferred]</b> Type of the \p identity value used Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <
        typename            InputIteratorRA,
        typename            OutputIteratorRA,
        typename            ScanOp,
        typename            Identity>
    __host__ __device__ __forceinline__
    static cudaError_t ExclusiveScan(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_in,                               ///< [in] Iterator pointing to scan input
        OutputIteratorRA    d_out,                              ///< [in] Iterator pointing to scan output
        ScanOp              scan_op,                            ///< [in] Binary scan operator
        Identity            identity,                           ///< [in] Identity element
        int                 num_items,                          ///< [in] Total number of items to scan
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch(d_temp_storage, temp_storage_bytes, d_in, d_out, scan_op, identity, num_items, stream, stream_synchronous);
    }


    //@}  end member group
    /******************************************************************//**
     * \name Inclusive scans
     *********************************************************************/
    //@{


    /**
     * \brief Computes a device-wide inclusive prefix sum.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the inclusive prefix sum of a device vector of \p int items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input and output
     * int *d_scan_input, *d_scan_output;
     * int num_items = ...
     * ...
     *
     * // Determine temporary device storage requirements for inclusive prefix sum
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, num_items);
     *
     * // Allocate temporary storage for inclusive prefix sum
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Run inclusive prefix sum
     * cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, num_items);
     *
     * \endcode
     *
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)
     * \tparam OutputIteratorRA     <b>[inferred]</b> Random-access iterator type for output (may be a simple pointer type)
     */
    template <
        typename            InputIteratorRA,
        typename            OutputIteratorRA>
    __host__ __device__ __forceinline__
    static cudaError_t InclusiveSum(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_in,                               ///< [in] Iterator pointing to scan input
        OutputIteratorRA    d_out,                              ///< [in] Iterator pointing to scan output
        int                 num_items,                          ///< [in] Total number of items to scan
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch(d_temp_storage, temp_storage_bytes, d_in, d_out, Sum(), NullType(), num_items, stream, stream_synchronous);
    }


    /**
     * \brief Computes a device-wide inclusive prefix scan using the specified binary \p scan_op functor.
     *
     * \par
     * Supports non-commutative scan operators.
     *
     * \devicestorage
     *
     * \cdp
     *
     * \iterator
     *
     * \par
     * The code snippet below illustrates the inclusive prefix scan of a device vector of \p int items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     * ...
     *
     * // Declare and initialize device pointers for input and output
     * int *d_scan_input, *d_scan_output;
     * int num_items = ...
     * ...
     *
     * // Determine temporary device storage requirements for inclusive prefix scan
     * void *d_temp_storage = NULL;
     * size_t temp_storage_bytes = 0;
     * cub::DeviceScan::InclusiveScan(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, cub::Max(), num_items);
     *
     * // Allocate temporary storage for inclusive prefix scan
     * cudaMalloc(&d_temp_storage, temp_storage_bytes);
     *
     * // Run inclusive prefix scan (max)
     * cub::DeviceScan::InclusiveScan(d_temp_storage, temp_storage_bytes, d_scan_input, d_scan_output, cub::Max(), num_items);
     *
     * \endcode
     *
     * \tparam InputIteratorRA      <b>[inferred]</b> Random-access iterator type for input (may be a simple pointer type)
     * \tparam OutputIteratorRA     <b>[inferred]</b> Random-access iterator type for output (may be a simple pointer type)
     * \tparam ScanOp               <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <
        typename            InputIteratorRA,
        typename            OutputIteratorRA,
        typename            ScanOp>
    __host__ __device__ __forceinline__
    static cudaError_t InclusiveScan(
        void                *d_temp_storage,                    ///< [in] %Device allocation of temporary storage.  When NULL, the required allocation size is returned in \p temp_storage_bytes and no work is done.
        size_t              &temp_storage_bytes,                ///< [in,out] Size in bytes of \p d_temp_storage allocation.
        InputIteratorRA     d_in,                               ///< [in] Iterator pointing to scan input
        OutputIteratorRA    d_out,                              ///< [in] Iterator pointing to scan output
        ScanOp              scan_op,                            ///< [in] Binary scan operator
        int                 num_items,                          ///< [in] Total number of items to scan
        cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                stream_synchronous  = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        return Dispatch(d_temp_storage, temp_storage_bytes, d_in, d_out, scan_op, NullType(), num_items, stream, stream_synchronous);
    }

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)


