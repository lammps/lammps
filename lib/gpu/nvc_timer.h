/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
------------------------------------------------------------------------- */

#ifndef NVC_TIMER_H
#define NVC_TIMER_H

#include "nvc_macros.h"

#define cudaEventDestroy(a)

/// Class for timing CUDA events
class NVCTimer {
 public:
  NVCTimer() : _total_time(0.0f), initialized(false) { }
  
  ~NVCTimer() {
    if (initialized) 
      { cudaEventDestroy(start_event); cudaEventDestroy(stop_event); }
  }

  inline void init() {
    if (initialized) {
      cudaEventDestroy(start_event); 
      cudaEventDestroy(stop_event);
    }
    initialized=true;
    CUDA_SAFE_CALL( cudaEventCreate(&start_event) );
    CUDA_SAFE_CALL( cudaEventCreate(&stop_event) );
  }
  
  /// Start timing
  inline void start() { cudaEventRecord(start_event,0); }
  
  /// Stop timing and store event time
  inline void stop() { cudaEventRecord(stop_event,0); }
  
  /// Set the time elapsed to zero (not the total_time)
  inline void zero() 
    { cudaEventRecord(start_event,0); cudaEventRecord(stop_event,0); } 
  
  /// Add time from previous start and stop to total
  /** Forces synchronization **/
  inline void add_to_total() { _total_time+=time(); }
  
  /// Return the time (ms) of last start to stop - Forces synchronization
  inline double time() { 
    float timer;
    cudaEventSynchronize(stop_event);
    CUDA_SAFE_CALL( cudaEventElapsedTime(&timer,start_event,stop_event) );
    return timer; 
  }
  
  /// Return the total time in ms
  inline double total_time() { return _total_time; }

  /// Return the total time in seconds
  inline double total_seconds() { return _total_time/1000.0; }

 private:
  cudaEvent_t start_event, stop_event;
  double _total_time;
  bool initialized;
};

#endif
