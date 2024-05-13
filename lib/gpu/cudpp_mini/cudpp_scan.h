// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 5289 $
// $Date: 2010-11-23 13:04:43 -0700 (Tue, 23 Nov 2010) $
// -------------------------------------------------------------
// This source code is distributed under the terms of license.txt
// in the root directory of this source distribution.
// -------------------------------------------------------------

/**
* @file
* cudpp_scan.h
*
* @brief Scan functionality header file - contains CUDPP interface (not public)
*/

#ifndef _CUDPP_SCAN_H_
#define _CUDPP_SCAN_H_

class CUDPPScanPlan;

extern "C"
void allocScanStorage(CUDPPScanPlan *plan);

extern "C"
void freeScanStorage(CUDPPScanPlan *plan);

extern "C"
void cudppScanDispatch(void                *d_out,
                       const void          *d_in,
                       size_t              numElements,
                       size_t              numRows,
                       const CUDPPScanPlan *plan);

#endif // _CUDPP_SCAN_H_
