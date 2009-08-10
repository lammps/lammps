/***************************************************************************
                                nvc_device.cu
                             -------------------
                               W. Michael Brown

  Utilities for dealing with cuda devices

 __________________________________________________________________________
    This file is part of the NVC Library
 __________________________________________________________________________

    begin                : Wed Jan 28 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "nvc_macros.h"
#include "nvc_device.h"

// Grabs the properties for all devices
NVCDevice::NVCDevice() {
  CUDA_SAFE_CALL(cudaGetDeviceCount(&_num_devices));
  for (int dev=0; dev<_num_devices; ++dev) {
    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));
    if (deviceProp.major == 9999 && deviceProp.minor == 9999)
      break;
    _properties.push_back(deviceProp);
  }
  _device=0;
}

// Set the CUDA device to the specified device number
void NVCDevice::set(int num) {
  if (_device==num)
    return;
  cudaThreadExit();
  CUDA_SAFE_CALL(cudaSetDevice(num));
  _device=num;
}

// List all devices along with all properties
void NVCDevice::print_all(ostream &out) {
  if (num_devices() == 0)
    printf("There is no device supporting CUDA\n");
  for (int i=0; i<num_devices(); ++i) {
    printf("\nDevice %d: \"%s\"\n", i, name(i).c_str());
    printf("  Revision number:                               %.1f\n", revision(i));
    printf("  Total amount of global memory:                 %.2f GB\n",
           gigabytes(i));
    #if CUDART_VERSION >= 2000
    printf("  Number of multiprocessors:                     %d\n",
           _properties[i].multiProcessorCount);
    printf("  Number of cores:                               %d\n",cores(i));
    #endif
    printf("  Total amount of constant memory:               %u bytes\n",
           _properties[i].totalConstMem); 
    printf("  Total amount of shared memory per block:       %u bytes\n",
           _properties[i].sharedMemPerBlock);
    printf("  Total number of registers available per block: %d\n",
           _properties[i].regsPerBlock);
    printf("  Warp size:                                     %d\n",
           _properties[i].warpSize);
    printf("  Maximum number of threads per block:           %d\n",
           _properties[i].maxThreadsPerBlock);
    printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
           _properties[i].maxThreadsDim[0],
           _properties[i].maxThreadsDim[1],
           _properties[i].maxThreadsDim[2]);
    printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
           _properties[i].maxGridSize[0],
           _properties[i].maxGridSize[1],
           _properties[i].maxGridSize[2]);
    printf("  Maximum memory pitch:                          %u bytes\n",
           _properties[i].memPitch);
    printf("  Texture alignment:                             %u bytes\n",
           _properties[i].textureAlignment);
    printf("  Clock rate:                                    %.2f GHz\n",
           clock_rate(i));
    #if CUDART_VERSION >= 2000
    printf("  Concurrent copy and execution:                 %s\n",
           _properties[i].deviceOverlap ? "Yes" : "No");
    #endif
  }
}

