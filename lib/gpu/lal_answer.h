/***************************************************************************
                                  answer.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for data management of forces, torques, energies, and virials

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_ANSWER_H
#define LAL_ANSWER_H

#include <math.h>
#include "mpi.h"

#ifdef USE_OPENCL

#include "geryon/ocl_timer.h"
#include "geryon/ocl_mat.h"
using namespace ucl_opencl;

#else

#include "geryon/nvd_timer.h"
#include "geryon/nvd_mat.h"
using namespace ucl_cudadr;

#endif

#include "lal_precision.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class Answer {
 public:
  Answer();
  ~Answer() { clear(); }

  /// Current number of local atoms stored
  inline int inum() const { return _inum; }
  /// Set number of local atoms for future copy operations
  inline void inum(const int n) { _inum=n; }
  
  /// Memory usage per atom in this class
  int bytes_per_atom() const; 

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param rot True if atom storage needs quaternions **/
  bool init(const int inum, const bool charge, const bool rot, UCL_Device &dev);
  
  /// Check if we have enough device storage and realloc if not
  inline void resize(const int inum, bool &success) {
    _inum=inum;
    if (inum>_max_local) {
      clear_resize();
      success = success && alloc(inum);
    }
  }
  
  /// If already initialized by another LAMMPS style, add fields as necessary
  /** \param rot True if atom storage needs quaternions **/
  bool add_fields(const bool charge, const bool rot);
  
  /// Free all memory on host and device needed to realloc for more atoms
  void clear_resize();

  /// Free all memory on host and device
  void clear();
 
  /// Return the total amount of host memory used by class in bytes
  double host_memory_usage() const;

  /// Add copy times to timers
  inline void acc_timers() {
    time_answer.add_to_total();
  }

  /// Add copy times to timers
  inline void zero_timers() {
    time_answer.zero();
  }

  /// Return the total time for host/device data transfer
  inline double transfer_time() {
    return time_answer.total_seconds();
  }
  
  /// Return the total time for data cast/pack
  inline double cast_time() { return _time_cast; }
  
  /// Return number of bytes used on device
  inline double gpu_bytes() { return _gpu_bytes; } 

  // -------------------------COPY FROM GPU -------------------------------

  /// Copy answers from device into read buffer asynchronously
  void copy_answers(const bool eflag, const bool vflag,
                    const bool ef_atom, const bool vf_atom);

  /// Copy answers from device into read buffer asynchronously
  void copy_answers(const bool eflag, const bool vflag,
                    const bool ef_atom, const bool vf_atom, int *ilist);
  
  /// Copy energy and virial data into LAMMPS memory
  double energy_virial(double *eatom, double **vatom, double *virial);

  /// Copy energy and virial data into LAMMPS memory
  double energy_virial(double *eatom, double **vatom, double *virial,
                       double &ecoul);

  /// Add forces and torques from the GPU into a LAMMPS pointer
  void get_answers(double **f, double **tor);

  inline double get_answers(double **f, double **tor, double *eatom, 
                            double **vatom, double *virial, double &ecoul) {
    double ta=MPI_Wtime();
    time_answer.sync_stop();
    _time_cpu_idle+=MPI_Wtime()-ta;
    double ts=MPI_Wtime();
    double evdw=energy_virial(eatom,vatom,virial,ecoul);
    get_answers(f,tor);
    _time_cast+=MPI_Wtime()-ts;
    return evdw;
  }
  
  /// Return the time the CPU was idle waiting for GPU
  inline double cpu_idle_time() { return _time_cpu_idle; }

  // ------------------------------ DATA ----------------------------------

  /// Force and possibly torque
  UCL_D_Vec<acctyp> dev_ans;
  /// Energy and virial per-atom storage
  UCL_D_Vec<acctyp> dev_engv;
  
  /// Force and possibly torque data on host
  UCL_H_Vec<acctyp> host_ans;
  /// Energy/virial data on host
  UCL_H_Vec<acctyp> host_engv;
  
  /// Device timers
  UCL_Timer time_answer;
  
  /// Geryon device
  UCL_Device *dev;

 private:
  bool alloc(const int inum);
  
  bool _allocated, _eflag, _vflag, _ef_atom, _vf_atom, _rot, _charge, _other;
  int _max_local, _inum, _e_fields, _ev_fields;
  int *_ilist;
  double _time_cast, _time_cpu_idle;
  
  double _gpu_bytes;
  
  bool _newton;
};

}

#endif
