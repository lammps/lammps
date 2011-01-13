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

#ifndef PAIR_GPU_ATOM_H
#define PAIR_GPU_ATOM_H

#include <math.h>
#include "mpi.h"

#ifdef USE_OPENCL

#include "geryon/ocl_device.h"
#include "geryon/ocl_timer.h"
#include "geryon/ocl_mat.h"
#include "geryon/ocl_kernel.h"
using namespace ucl_opencl;

#else

#include "cudpp.h"
#include "geryon/nvd_device.h"
#include "geryon/nvd_timer.h"
#include "geryon/nvd_mat.h"
#include "geryon/nvd_kernel.h"
using namespace ucl_cudadr;

#endif

#ifndef int2
struct int2 { int x; int y; };
#endif

#include "pair_gpu_precision.h"

template <class numtyp, class acctyp>
class PairGPUAtom {
 public:
  PairGPUAtom();
  ~PairGPUAtom() { clear(); }

  /// Maximum number of atoms that can be stored with current allocation
  inline int max_atoms() const { return _max_atoms; }
  /// Current number of local+ghost atoms stored
  inline int nall() const { return _nall; }
  /// Current number of local atoms stored
  inline int inum() const { return _inum; }

  /// Set number of local+ghost atoms for future copy operations
  inline void nall(const int n) { _nall=n; }
  /// Set number of local atoms for future copy operations
  inline void inum(const int n) { _inum=n; }
  
  /// Memory usage per atom in this class
  int bytes_per_atom() const; 

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param rot True if atom storage needs quaternions
    * \param gpu_nbor True if neighboring will be performed on device **/
  bool init(const int inum, const int nall, const bool charge, const bool rot, 
            UCL_Device &dev, const bool gpu_nbor=false, const bool bonds=false);
  
  /// Check if we have enough device storage and realloc if not
  inline bool resize(const int inum, const int nall, bool &success) {
    _inum=inum;
    _nall=nall;
    if (inum>_max_local || nall>_max_atoms) {
      clear_resize();
      success = success && alloc(inum,nall);
      return true;
    }
    return false;
  }

  /// Only free matrices of length inum or nall for resizing
  void clear_resize();
  
  /// Free all memory on host and device
  void clear();
 
  /// Return the total amount of host memory used by class in bytes
  double host_memory_usage() const;

  /// Sort arrays for neighbor list calculation on device
  void sort_neighbor(const int num_atoms);
  
  /// Add copy times to timers
  inline void acc_timers() {
    time_pos.add_to_total();
    time_answer.add_to_total();
    if (_other)
      time_other.add_to_total();
  }

  /// Add copy times to timers
  inline void zero_timers() {
    time_pos.zero();
    time_answer.zero();
    if (_other)
      time_other.zero();
  }

  /// Return the total time for host/device data transfer
  inline double transfer_time() {
    double total=time_pos.total_seconds()+time_answer.total_seconds();
    if (_other) total+=time_other.total_seconds();
    return total;
  }
  
  /// Return the total time for data cast/pack
  inline double cast_time() { return _time_cast; }

  /// Pack LAMMPS atom type constants into matrix and copy to device
  template <class dev_typ, class t1>
  inline void type_pack1(const int n, const int m_size,
			 UCL_D_Vec<dev_typ> &dev_v, UCL_H_Vec<numtyp> &buffer,
			 t1 **one) {
    int ii=0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        buffer[ii]=static_cast<numtyp>(one[i][j]);
        ii++;
      }
      ii+=m_size-n;
    }
    UCL_H_Vec<dev_typ> view;
    view.view((dev_typ*)buffer.begin(),m_size*m_size,*dev);
    ucl_copy(dev_v,view,false);
  }

  /// Pack LAMMPS atom type constants into 2 vectors and copy to device
  template <class dev_typ, class t1, class t2>
  inline void type_pack2(const int n, const int m_size,
			 UCL_D_Vec<dev_typ> &dev_v, UCL_H_Vec<numtyp> &buffer,
			 t1 **one, t2 **two) {
    int ii=0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        buffer[ii*2]=static_cast<numtyp>(one[i][j]);
        buffer[ii*2+1]=static_cast<numtyp>(two[i][j]);
        ii++;
      }
      ii+=m_size-n;
    }
    UCL_H_Vec<dev_typ> view;
    view.view((dev_typ*)buffer.begin(),m_size*m_size,*dev);
    ucl_copy(dev_v,view,false);
  }

  /// Pack LAMMPS atom type constants (3) into 4 vectors and copy to device
  template <class dev_typ, class t1, class t2, class t3>
  inline void type_pack4(const int n, const int m_size,
			 UCL_D_Vec<dev_typ> &dev_v, UCL_H_Vec<numtyp> &buffer,
			 t1 **one, t2 **two, t3 **three) {
    int ii=0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        buffer[ii*4]=static_cast<numtyp>(one[i][j]);
        buffer[ii*4+1]=static_cast<numtyp>(two[i][j]);
        buffer[ii*4+2]=static_cast<numtyp>(three[i][j]);
        ii++;
      }
      ii+=m_size-n;
    }
    UCL_H_Vec<dev_typ> view;
    view.view((dev_typ*)buffer.begin(),m_size*m_size,*dev);
    ucl_copy(dev_v,view,false);
  }
  
  /// Pack LAMMPS atom type constants (4) into 4 vectors and copy to device
  template <class dev_typ, class t1, class t2, class t3, class t4>
  inline void type_pack4(const int n, const int m_size,
			 UCL_D_Vec<dev_typ> &dev_v, UCL_H_Vec<numtyp> &buffer,
			 t1 **one, t2 **two, t3 **three, t4 **four) {
    int ii=0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        buffer[ii*4]=static_cast<numtyp>(one[i][j]);
        buffer[ii*4+1]=static_cast<numtyp>(two[i][j]);
        buffer[ii*4+2]=static_cast<numtyp>(three[i][j]);
        buffer[ii*4+3]=static_cast<numtyp>(four[i][j]);
        ii++;
      }
      ii+=m_size-n;
    }
    UCL_H_Vec<dev_typ> view;
    view.view((dev_typ*)buffer.begin(),m_size*m_size,*dev);
    ucl_copy(dev_v,view,false);
  }

  /// Pack LAMMPS atom "self" type constants into 2 vectors and copy to device
  template <class dev_typ, class t1, class t2>
  inline void self_pack2(const int n, UCL_D_Vec<dev_typ> &dev_v, 
                         UCL_H_Vec<numtyp> &buffer, t1 **one, t2 **two) {
    for (int i=0; i<n; i++) {
      buffer[i*2]=static_cast<numtyp>(one[i][i]);
      buffer[i*2+1]=static_cast<numtyp>(two[i][i]);
    }
    UCL_H_Vec<dev_typ> view;
    view.view((dev_typ*)buffer.begin(),n,*dev);
    ucl_copy(dev_v,view,false);
  }

  // -------------------------COPY TO GPU ----------------------------------

  /// Cast positions and types to write buffer
  inline void cast_x_data(double **host_ptr, const int *host_type) {
    double t=MPI_Wtime();
    #ifdef GPU_CAST
    memcpy(host_x_cast.begin(),host_ptr[0],_nall*3*sizeof(double));
    memcpy(host_type_cast.begin(),host_type,_nall*sizeof(int));
    #else
    numtyp *_write_loc=host_x.begin();
    for (int i=0; i<_nall; i++) {
      *_write_loc=host_ptr[i][0];
      _write_loc++;
      *_write_loc=host_ptr[i][1];
      _write_loc++;
      *_write_loc=host_ptr[i][2];
      _write_loc++;
      *_write_loc=host_type[i];
      _write_loc++;
    }
    #endif
    _time_cast+=MPI_Wtime()-t;
  }      

  /// Copy positions and types to device asynchronously
  /** Copies nall() elements **/
  inline void add_x_data(double **host_ptr, int *host_type) { 
    time_pos.start();
    #ifdef GPU_CAST
    ucl_copy(dev_x_cast,host_x_cast,_nall*3,true);
    ucl_copy(dev_type_cast,host_type_cast,_nall,true);
    int block_size=64;
    int GX=static_cast<int>(ceil(static_cast<double>(_nall)/block_size));
    k_cast_x.set_size(GX,block_size);
    k_cast_x.run(&dev_x.begin(), &dev_x_cast.begin(), &dev_type_cast.begin(), 
                 &_nall);
    #else
    ucl_copy(dev_x,host_x,_nall*4,true);
    #endif
    time_pos.stop();
  }

  /// Calls cast_x_data and add_x_data and times the routines
  inline void cast_copy_x(double **host_ptr, int *host_type) {
    cast_x_data(host_ptr,host_type);
    add_x_data(host_ptr,host_type);
  }

  /// Cast charges to write buffer
  template<class cpytyp>
  inline void cast_q_data(cpytyp *host_ptr) {
    double t=MPI_Wtime();
    if (dev->device_type()==UCL_CPU) {
      if (sizeof(numtyp)==sizeof(double)) {
        host_q.view((numtyp*)host_ptr,_nall,*dev);
        dev_q.view(host_q);
      } else
        for (int i=0; i<_nall; i++) host_q[i]=host_ptr[i];
    } else {
      if (sizeof(numtyp)==sizeof(double))
        memcpy(host_q.begin(),host_ptr,_nall*sizeof(numtyp));
      else
        for (int i=0; i<_nall; i++) host_q[i]=host_ptr[i];
    }
    _time_cast+=MPI_Wtime()-t;
  }

  /// Copy charges to device asynchronously
  inline void add_q_data() {
    ucl_copy(dev_q,host_q,_nall,true);
  }

  /// Cast quaternions to write buffer
  template<class cpytyp>
  inline void cast_quat_data(cpytyp *host_ptr) {
    double t=MPI_Wtime();
    if (dev->device_type()==UCL_CPU) {
      if (sizeof(numtyp)==sizeof(double)) {
        host_quat.view((numtyp*)host_ptr,_nall*4,*dev);
        dev_quat.view(host_quat);
      } else
        for (int i=0; i<_nall*4; i++) host_quat[i]=host_ptr[i];
    } else {
      if (sizeof(numtyp)==sizeof(double))
        memcpy(host_quat.begin(),host_ptr,_nall*4*sizeof(numtyp));
      else
        for (int i=0; i<_nall*4; i++) host_quat[i]=host_ptr[i];
    }
    _time_cast+=MPI_Wtime()-t;
  }

  /// Copy quaternions to device
  /** Copies nall()*4 elements **/
  inline void add_quat_data() {
    ucl_copy(dev_quat,host_quat,_nall*4,true);
  }

  /// Copy data other than pos and data to device
  inline void add_other_data() {
    time_other.start();
    if (_charge)
      add_q_data();
    if (_rot)
      add_quat_data();
    time_other.stop();
  }
  
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

  // ------------------------------ DATA ----------------------------------

  /// Atom coordinates and types ([0] is x, [1] is y, [2] is z, [3] is type
  UCL_D_Vec<numtyp> dev_x;
  /// Charges
  UCL_D_Vec<numtyp> dev_q;
  /// Quaterions
  UCL_D_Vec<numtyp> dev_quat;
  /// Force and possibly torque
  UCL_D_Vec<acctyp> dev_ans;
  /// Energy and virial per-atom storage
  UCL_D_Vec<acctyp> dev_engv;
  
  #ifdef GPU_CAST
  UCL_D_Vec<double> dev_x_cast;
  UCL_D_Vec<int> dev_type_cast;
  UCL_H_Vec<double> host_x_cast;
  UCL_H_Vec<int> host_type_cast;
  #endif

  /// Buffer for moving positions to device
  UCL_H_Vec<numtyp> host_x;
  /// Buffer for moving charge data to GPU
  UCL_H_Vec<numtyp> host_q;
  /// Buffer for moving quat data to GPU
  UCL_H_Vec<numtyp> host_quat;
  /// Force and possibly torque data on host
  UCL_H_Vec<acctyp> host_ans;
  /// Energy/virial data on host
  UCL_H_Vec<acctyp> host_engv;
  
  /// Cell list identifiers for device nbor builds
  UCL_D_Vec<unsigned> dev_cell_id;
  /// Cell list identifiers for device nbor builds
  UCL_D_Vec<int> dev_particle_id;
  /// Atom tag information for device nbor builds
  UCL_D_Vec<int> dev_tag;

  /// Device timers
  UCL_Timer time_pos, time_other, time_answer;
  
  /// Geryon device
  UCL_Device *dev;

 private:
  #ifdef GPU_CAST
  UCL_Program *atom_program;
  UCL_Kernel k_cast_x;
  void compile_kernels(UCL_Device &dev);
  #endif

  bool _compiled;

  bool alloc(const int inum, const int nall);
  
  bool _allocated, _eflag, _vflag, _ef_atom, _vf_atom, _rot, _charge, _other;
  int _max_local, _max_atoms, _nall, _inum, _e_fields, _ev_fields;
  bool _gpu_nbor, _bonds;
  int *_ilist;
  double _time_cast;
  
  double _gpu_bytes;
  
  bool _newton;

  #ifndef USE_OPENCL
  CUDPPConfiguration sort_config;
  CUDPPHandle sort_plan;
  #endif
};

#endif

