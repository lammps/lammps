// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision$
// $Date$
// -------------------------------------------------------------
// This source code is distributed under the terms of license.txt
// in the root directory of this source distribution.
// -------------------------------------------------------------
#include "cudpp_maximal_launch.h"
#include <algorithm>

// computes next highest multiple of f from x
inline size_t multiple(size_t x, size_t f)
{
    return ((x + (f-1)) / f);
}


// MS Excel-style CEIL() function
// Rounds x up to nearest multiple of f
inline size_t ceiling(size_t x, size_t f)
{
    return multiple(x, f) * f;
}

extern "C"
size_t maxBlocks(cudaFuncAttributes &attribs,
                 cudaDeviceProp &devprop,
                 size_t bytesDynamicSharedMem,
                 size_t threadsPerBlock)
{

    // Determine the maximum number of CTAs that can be run simultaneously for each kernel
    // This is equivalent to the calculation done in the CUDA Occupancy Calculator spreadsheet
    const unsigned int regAllocationUnit = (devprop.major < 2 && devprop.minor < 2) ? 256 : 512; // in registers
    const unsigned int warpAllocationMultiple = 2;
    const unsigned int smemAllocationUnit = 512;                                                 // in bytes
    const unsigned int maxThreadsPerSM = (devprop.major < 2 && devprop.minor < 2) ? 768 : 1024;  // sm_12 GPUs increase threads/SM to 1024
    const unsigned int maxBlocksPerSM = 8;

    // Number of warps (round up to nearest whole multiple of warp size)
    size_t numWarps = multiple(threadsPerBlock, devprop.warpSize);
    // Round up to warp allocation multiple
    numWarps = ceiling(numWarps, warpAllocationMultiple);

    // Number of regs is regs per thread times number of warps times warp size
    size_t regsPerCTA = attribs.numRegs * devprop.warpSize * numWarps;
    // Round up to multiple of register allocation unit size
    regsPerCTA = ceiling(regsPerCTA, regAllocationUnit);

    size_t smemBytes  = attribs.sharedSizeBytes + bytesDynamicSharedMem;
    size_t smemPerCTA = ceiling(smemBytes, smemAllocationUnit);

    size_t ctaLimitRegs    = regsPerCTA > 0 ? devprop.regsPerBlock      / regsPerCTA : maxBlocksPerSM;
    size_t ctaLimitSMem    = smemPerCTA > 0 ? devprop.sharedMemPerBlock / smemPerCTA : maxBlocksPerSM;
    size_t ctaLimitThreads =                  maxThreadsPerSM           / threadsPerBlock;

    return devprop.multiProcessorCount * std::min(ctaLimitRegs, std::min(ctaLimitSMem, std::min(ctaLimitThreads, (size_t)maxBlocksPerSM)));
}

extern "C"
size_t maxBlocksFromPointer(void*  kernel,
                            size_t bytesDynamicSharedMem,
                            size_t threadsPerBlock)
{
    cudaDeviceProp devprop;
    int deviceID = -1;
    cudaError_t err = cudaGetDevice(&deviceID);
    if (err == cudaSuccess)
    {
        err = cudaGetDeviceProperties(&devprop, deviceID);
        if (err != cudaSuccess)
            return (size_t)-1;

        cudaFuncAttributes attr;
        err = cudaFuncGetAttributes(&attr, (const char*)kernel);
        if (err != cudaSuccess)
            return (size_t)-1;

        return maxBlocks(attr, devprop, bytesDynamicSharedMem, threadsPerBlock);
    }

    return (size_t)-1;
}
