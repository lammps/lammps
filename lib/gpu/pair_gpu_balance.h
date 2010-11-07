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
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifndef PAIR_GPU_BALANCE_H
#define PAIR_GPU_BALANCE_H

#include "pair_gpu_device.h"
#include <math.h>

#define _HD_BALANCE_EVERY 25
#define _HD_BALANCE_WEIGHT 0.5
#define _HD_BALANCE_GAP 1.05

/// Host/device load balancer
template<class numtyp, class acctyp>
class PairGPUBalance {
 public:
  inline PairGPUBalance() : _init_done(false), _measure_this_step(false) {}
  inline ~PairGPUBalance() { clear(); }

  /// Clear any old data and setup for new LAMMPS run
  inline void init(PairGPUDevice<numtyp, acctyp> *gpu, const double split);

  /// Clear all host and device data
  inline void clear() {
    if (_init_done) {
      _device_time.clear();
      _measure_this_step=false;
      _init_done=false;
    }
  }

  /// Get a count of the number of particles host will handle for initial alloc
  inline int first_host_count(const int nlocal,const bool gpu_nbor,
                              const double gpu_split) const {
    int host_nlocal=0;
    if (gpu_nbor && gpu_split!=1.0) {
      if (gpu_split>0)
        host_nlocal=static_cast<int>(ceil((1.0-gpu_split)*nlocal));
      else
        host_nlocal=static_cast<int>(ceil(0.1*nlocal));
    }
    return host_nlocal;
  }

  /// Return the number of particles the device will handle this timestep
  inline int get_gpu_count(const int timestep, const int ago,
                           const int inum_full);

  /// Return the average fraction of particles handled by device on all procs
  inline double all_avg_split() {
    if (_load_balance) {
      double _all_avg_split=0.0;
      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
      MPI_Reduce(&_avg_split,&_all_avg_split,1,MPI_DOUBLE,MPI_SUM,0,
                 MPI_COMM_WORLD);
      _all_avg_split/=nprocs;
      return _all_avg_split/_avg_count;
    } else
      return _actual_split;
  }

  /// If CPU neighboring, allow the device fraction to increase on 2nd timestep
  inline int ago_first(int ago) const
    { if (_avg_count==1 && _actual_split<_desired_split) ago=0; return ago; }

  /// Start the timer for asynchronous device execution
  inline void start_timer() {
    if (_measure_this_step) {
      _device->gpu->sync();
      MPI_Barrier(_device->gpu_comm);
      _device_time.start();
      _device->gpu->sync();
      MPI_Barrier(_device->gpu_comm);
      _device->start_host_timer();
    }
  }

  /// Stop the timer for asynchronous device execution
  inline void stop_timer() { if (_measure_this_step) { _device_time.stop(); } }

  /// Calculate the new host/device split based on the cpu and device times
  /** \note Only does calculation every _HD_BALANCE_EVERY timesteps 
            (and first 10) **/
  inline void balance(const double cpu_time, const bool gpu_nbor);

  /// Calls balance() and then get_gpu_count()
  inline int balance(const int timestep, const int ago, const int inum_full,
                     const double cpu_time, const bool gpu_nbor) {
    balance(cpu_time,gpu_nbor);
    return get_gpu_count(timestep,ago,inum_full);
  }
  
 private:
  PairGPUDevice<numtyp,acctyp> *_device;
  UCL_Timer _device_time;
  bool _init_done;
  
  bool _load_balance;
  double _actual_split, _avg_split, _desired_split, _max_split;
  int _avg_count;

  bool _measure_this_step;
  int _inum, _inum_full;
};

#define PairGPUBalanceT PairGPUBalance<numtyp,acctyp>

template <class numtyp, class acctyp>
void PairGPUBalanceT::init(PairGPUDevice<numtyp, acctyp> *gpu,
			   const double split) {
  clear();
  _init_done=true;
  
  _device=gpu;
  _device_time.init(*gpu->gpu);
  
  if (split<0.0) {
    _load_balance=true;
    _desired_split=0.9;
  } else {
    _load_balance=false;
    _desired_split=split;
  }
  _actual_split=_desired_split;
  _avg_split=0.0;
  _avg_count=0;
}

template <class numtyp, class acctyp>
int PairGPUBalanceT::get_gpu_count(const int timestep, const int ago,
			           const int inum_full) {
  _measure_this_step=false;
  if (_load_balance) {
    if (_avg_count<11 || timestep%_HD_BALANCE_EVERY==0) {
      _measure_this_step=true;
      _inum_full=inum_full;
    }
    if (ago==0) {
      _actual_split=_desired_split;
      _max_split=_desired_split;
    }
  }
  _inum=static_cast<int>(floor(_actual_split*inum_full));
  if (_inum==0) _inum++;
  return _inum;
}
    
template <class numtyp, class acctyp>
void PairGPUBalanceT::balance(const double cpu_time, const bool gpu_nbor) {
  if (_measure_this_step) {
    if (_inum_full==_inum) {
      _desired_split=1.0;
      return;
    }

    _measure_this_step=false;
    double gpu_time=_device_time.seconds();

    double cpu_gpu_time[3], max_times[3];
    cpu_gpu_time[0]=cpu_time/(_inum_full-_inum);
    cpu_gpu_time[1]=gpu_time/_inum;
    cpu_gpu_time[2]=(_device->host_time()-cpu_time)/_inum_full;

    MPI_Allreduce(cpu_gpu_time,max_times,3,MPI_DOUBLE,MPI_MAX,
                  _device->gpu_comm);
    double split=(max_times[0]+max_times[2])/(max_times[0]+max_times[1]);
    split*=_HD_BALANCE_GAP;

    if (split>1.0)
      split=1.0;
    if (_avg_count<10)
      _desired_split=(_desired_split*_avg_count+split)/(_avg_count+1);
    else
      _desired_split=_desired_split*(1.0-_HD_BALANCE_WEIGHT)+
                     _HD_BALANCE_WEIGHT*split;

    if (!gpu_nbor) {
      if (_desired_split<_max_split)
        _actual_split=_desired_split;
      else
        _actual_split=_max_split;
    }
  }
  _avg_split+=_desired_split;
  _avg_count++;
}

#endif

