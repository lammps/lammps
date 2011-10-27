/***************************************************************************
                                  balance.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for host-device load balancing

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_BALANCE_H
#define LAL_BALANCE_H

#include "lal_device.h"
#include <math.h>

#define _HD_BALANCE_EVERY 25
#define _HD_BALANCE_WEIGHT 0.5
#define _HD_BALANCE_GAP 1.10

namespace LAMMPS_AL {

/// Host/device load balancer
template<class numtyp, class acctyp>
class Balance {
 public:
  inline Balance() : _init_done(false), _measure_this_step(false) {}
  inline ~Balance() { clear(); }

  /// Clear any old data and setup for new LAMMPS run
  inline void init(Device<numtyp, acctyp> *gpu, const int gpu_nbor,
                   const double split);

  /// Clear all host and device data
  inline void clear() {
    if (_init_done) {
      _device_time.clear();
      _measure_this_step=false;
      _init_done=false;
    }
  }
  
  /// Return the timestep since initialization
  inline int timestep() { return _timestep; }

  /// Get a count of the number of particles host will handle for initial alloc
  inline int first_host_count(const int nlocal, const double gpu_split,
                              const int gpu_nbor) const {
    int host_nlocal=0;
    if (gpu_nbor>0 && gpu_split!=1.0) {
      if (gpu_split>0)
        host_nlocal=static_cast<int>(ceil((1.0-gpu_split)*nlocal));
      else
        host_nlocal=static_cast<int>(ceil(0.05*nlocal));
    }
    return host_nlocal;
  }

  /// Return the number of particles the device will handle this timestep
  inline int get_gpu_count(const int ago, const int inum_full);

  /// Return the average fraction of particles handled by device on all procs
  inline double all_avg_split() {
    if (_load_balance) {
      double _all_avg_split=0.0;
      MPI_Reduce(&_avg_split,&_all_avg_split,1,MPI_DOUBLE,MPI_SUM,0,
                 _device->replica());
      _all_avg_split/=_device->replica_size();
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
      _device->gpu_barrier();
      _device->start_host_timer();
      _device_time.start();
      _device->gpu->sync();
      _device->gpu_barrier();
    }
  }

  /// Stop the timer for asynchronous device execution
  inline void stop_timer() { if (_measure_this_step) { _device_time.stop(); } }

  /// Calculate the new host/device split based on the cpu and device times
  /** \note Only does calculation every _HD_BALANCE_EVERY timesteps 
            (and first 10) **/
  inline void balance(const double cpu_time);

  /// Calls balance() and then get_gpu_count()
  inline int balance(const int ago,const int inum_full,const double cpu_time) {
    balance(cpu_time);
    return get_gpu_count(ago,inum_full);
  }
  
 private:
  Device<numtyp,acctyp> *_device;
  UCL_Timer _device_time;
  bool _init_done;
  int _gpu_nbor;
  
  bool _load_balance;
  double _actual_split, _avg_split, _desired_split, _max_split;
  int _avg_count;

  bool _measure_this_step;
  int _inum, _inum_full, _timestep;
};

#define BalanceT Balance<numtyp,acctyp>

template <class numtyp, class acctyp>
void BalanceT::init(Device<numtyp, acctyp> *gpu, 
                           const int gpu_nbor, const double split) {
  clear();
  _gpu_nbor=gpu_nbor;
  _init_done=true;
  
  _device=gpu;
  _device_time.init(*gpu->gpu);
  
  if (split<0.0) {
    _load_balance=true;
    _desired_split=0.90;
  } else {
    _load_balance=false;
    _desired_split=split;
  }
  _actual_split=_desired_split;
  _avg_split=0.0;
  _avg_count=0;
  _timestep=0;
}

template <class numtyp, class acctyp>
int BalanceT::get_gpu_count(const int ago, const int inum_full) {
  _measure_this_step=false;
  if (_load_balance) {
    if (_avg_count<11 || _timestep%_HD_BALANCE_EVERY==0) {
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
  _timestep++;
  return _inum;
}
    
template <class numtyp, class acctyp>
void BalanceT::balance(const double cpu_time) {
  if (_measure_this_step) {
    _measure_this_step=false;
    double gpu_time=_device_time.seconds();

    double max_gpu_time;
    MPI_Allreduce(&gpu_time,&max_gpu_time,1,MPI_DOUBLE,MPI_MAX,
                  _device->gpu_comm());

    if (_inum_full==_inum) {
      _desired_split=1.0;
      return;
    }

    double cpu_time_per_atom=cpu_time/(_inum_full-_inum);
    double cpu_other_time=_device->host_time()-cpu_time;
    int host_inum=static_cast<int>((max_gpu_time-cpu_other_time)/
                                   cpu_time_per_atom);

    double split=static_cast<double>(_inum_full-host_inum)/_inum_full;
    _desired_split=split*_HD_BALANCE_GAP;
    if (_desired_split>1.0)
      _desired_split=1.0;
    if (_desired_split<0.0)
      _desired_split=0.0;

    if (_gpu_nbor==0) {
      if (_desired_split<_max_split)
        _actual_split=_desired_split;
      else
        _actual_split=_max_split;
    }
  }
  _avg_split+=_desired_split;
  _avg_count++;
}

}

#endif
