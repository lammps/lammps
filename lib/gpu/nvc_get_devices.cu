/***************************************************************************
                              nvc_get_devices.h
                             -------------------
                               W. Michael Brown

  List properties of cuda devices

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

#include "nvc_device.h"

int main(int argc, char** argv) {
  NVCDevice gpu;
  gpu.print_all(cout);
  return 0;
}

