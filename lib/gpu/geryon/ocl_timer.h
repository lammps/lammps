/***************************************************************************
                                 ocl_timer.h
                             -------------------
                               W. Michael Brown

  Class for timing OpenCL routines

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Jan Fri 22 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef OCL_TIMER_H
#define OCL_TIMER_H

#include "ocl_macros.h"
#include "ocl_device.h"

#ifdef CL_VERSION_1_2
#define UCL_OCL_MARKER(cq,event) clEnqueueMarkerWithWaitList(cq,0,nullptr,event)
#else
#define UCL_OCL_MARKER clEnqueueMarker
#endif

namespace ucl_opencl {

/// Class for timing OpenCL events
class UCL_Timer {
 public:
  inline UCL_Timer() : start_event(nullptr), stop_event(nullptr), _total_time(0.0f),
                       _initialized(false), has_measured_time(false) { }
  inline UCL_Timer(UCL_Device &dev) : start_event(nullptr), stop_event(nullptr), _total_time(0.0f),
                                      _initialized(false), has_measured_time(false)
    { init(dev); }

  inline ~UCL_Timer() { clear(); }

  /// Clear any data associated with timer
  /** \note init() must be called to reuse timer after a clear() **/
  inline void clear() {
    if (_initialized) {
      CL_DESTRUCT_CALL(clReleaseCommandQueue(_cq));
      _initialized=false;
      _total_time=0.0;
    }
    has_measured_time = false;
  }

  /// Initialize default command queue for timing
  inline void init(UCL_Device &dev) { init(dev,dev.cq()); }

  /// Initialize command queue for timing
  inline void init(UCL_Device &dev, command_queue &cq) {
    clear();
    _cq=cq;
    clRetainCommandQueue(_cq);
    _initialized=true;
    has_measured_time = false;
  }

  /// Start timing on default command queue
  inline void start() {
    UCL_OCL_MARKER(_cq,&start_event);
    has_measured_time = false;
  }

  /// Stop timing on default command queue
  inline void stop() {
    UCL_OCL_MARKER(_cq,&stop_event);
    has_measured_time = true;
  }

  /// Block until the start event has been reached on device
  inline void sync_start() {
    CL_SAFE_CALL(clWaitForEvents(1,&start_event));
    has_measured_time = false;
  }

  /// Block until the stop event has been reached on device
  inline void sync_stop() {
    CL_SAFE_CALL(clWaitForEvents(1,&stop_event));
    has_measured_time = true;
  }

  /// Set the time elapsed to zero (not the total_time)
  inline void zero() {
    has_measured_time = false;
    UCL_OCL_MARKER(_cq,&start_event);
    UCL_OCL_MARKER(_cq,&stop_event);
  }

  /// Set the total time to zero
  inline void zero_total() { _total_time=0.0; }

  /// Add time from previous start and stop to total
  /** Forces synchronization **/
  inline double add_to_total()
    { double t=time(); _total_time+=t; return t/1000.0; }

  /// Add a user specified time to the total (ms)
  inline void add_time_to_total(const double t) { _total_time+=t; }

  /// Return the time (ms) of last start to stop - Forces synchronization
  inline double time() {
    if(!has_measured_time) return 0.0;
    cl_ulong tstart,tend;
    CL_SAFE_CALL(clWaitForEvents(1,&stop_event));
    CL_SAFE_CALL(clGetEventProfilingInfo(stop_event,
                                         CL_PROFILING_COMMAND_START,
                                         sizeof(cl_ulong), &tend, nullptr));
    CL_SAFE_CALL(clGetEventProfilingInfo(start_event,
                                         CL_PROFILING_COMMAND_END,
                                         sizeof(cl_ulong), &tstart, nullptr));
    clReleaseEvent(start_event);
    clReleaseEvent(stop_event);
    has_measured_time = false;
    return (tend-tstart)*1e-6;
  }

  /// Return the time (s) of last start to stop - Forces synchronization
  inline double seconds() { return time()*1e-3; }

  /// Return the total time in ms
  inline double total_time() { return _total_time; }

  /// Return the total time in seconds
  inline double total_seconds() { return _total_time*1e-3; }

 private:
  cl_event start_event, stop_event;
  cl_command_queue _cq;
  double _total_time;
  bool _initialized;
  bool has_measured_time;
};

} // namespace

#endif
