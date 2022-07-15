/***************************************************************************
                              nvc_get_devices.h
                             -------------------
                               W. Michael Brown

  List properties of cuda devices

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Wed Jan 28 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifdef UCL_OPENCL
#include "ocl_device.h"
using namespace ucl_opencl;
#endif

#ifdef UCL_CUDADR
#include "nvd_device.h"
using namespace ucl_cudadr;
#endif

#ifdef UCL_CUDART
#include "nvc_device.h"
using namespace ucl_cudart;
#endif

#ifdef UCL_HIP
#include "hip_device.h"
using namespace ucl_hip;
#endif

int main(int /*argc*/, char** /*argv*/) {
  UCL_Device cop;
  std::cout << "Found " << cop.num_platforms() << " platform(s).\n";
  if (cop.num_platforms()>0)
    cop.print_all(std::cout);
  return 0;
}

