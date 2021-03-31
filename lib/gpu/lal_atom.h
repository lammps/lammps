/***************************************************************************
                                   atom.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for particle data management

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef PAIR_GPU_ATOM_H
#define PAIR_GPU_ATOM_H

#include <cmath>
#include "mpi.h"

#if defined(USE_OPENCL)
#include "geryon/ocl_timer.h"
#include "geryon/ocl_mat.h"
#include "geryon/ocl_kernel.h"
using namespace ucl_opencl;
#ifndef LAL_NO_OCL_EV_JIT
#define LAL_OCL_EV_JIT
#endif
#elif defined(USE_CUDART)
#include "geryon/nvc_timer.h"
#include "geryon/nvc_mat.h"
#include "geryon/nvc_kernel.h"
using namespace ucl_cudart;
#elif defined(USE_HIP)
#include "geryon/hip_timer.h"
#include "geryon/hip_mat.h"
#include "geryon/hip_kernel.h"
using namespace ucl_hip;
#else
#include "geryon/nvd_timer.h"
#include "geryon/nvd_mat.h"
#include "geryon/nvd_kernel.h"
using namespace ucl_cudadr;
#endif

#ifdef USE_CUDPP
#include "cudpp.h"
#endif

#include "lal_precision.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class Atom {
 public:
  Atom();
  ~Atom() { clear(); }

  /// Maximum number of atoms that can be stored with current allocation
  inline int max_atoms() const { return _max_atoms; }
  /// Current number of local+ghost atoms stored
  inline int nall() const { return _nall; }

  /// Set number of local+ghost atoms for future copy operations
  inline void nall(const int n) { _nall=n; }

  /// Memory usage per atom in this class
  int bytes_per_atom() const;

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param rot True if atom storage needs quaternions
    * \param gpu_nbor 0 if neighboring will be performed on host
    *        gpu_nbor 1 if neighboring will be performed on device
    *        gpu_nbor 2 if binning on host and neighboring on device **/
  bool init(const int nall, const bool charge, const bool rot,
            UCL_Device &dev, const int gpu_nbor=0, const bool bonds=false,
            const bool vel=false);

  /// Check if we have enough device storage and realloc if not
  /** Returns true if resized with any call during this timestep **/
  inline bool resize(const int nall, bool &success) {
    _nall=nall;
    if (nall>_max_atoms) {
      clear_resize();
      success = success && alloc(nall);
      _resized=true;
    }
    return _resized;
  }

  /// If already initialized by another LAMMPS style, add fields as necessary
  /** \param rot True if atom storage needs quaternions
    * \param gpu_nbor 0 if neighboring will be performed on host
    *        gpu_nbor 1 if neighboring will be performed on device
    *        gpu_nbor 2 if binning on host and neighboring on device **/
  bool add_fields(const bool charge, const bool rot, const int gpu_nbor,
                  const bool bonds, const bool vel=false);

  /// Returns true if GPU is using charges
  bool charge() { return _charge; }

  /// Returns true if GPU is using quaternions
  bool quaternion() { return _rot; }

  /// Returns true if GPU is using velocities
  bool velocity() { return _vel; }

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
    if (_charge)
      time_q.add_to_total();
    if (_rot)
      time_quat.add_to_total();
    if (_vel)
      time_vel.add_to_total();
  }

  /// Add copy times to timers
  inline void zero_timers() {
    time_pos.zero();
    if (_charge)
      time_q.zero();
    if (_rot)
      time_quat.zero();
    if (_vel)
      time_vel.zero();
  }

  /// Return the total time for host/device data transfer
  /** Zeros the total so that the atom times are only included once **/
  inline double transfer_time() {
    double total=time_pos.total_seconds();
    time_pos.zero_total();
    if (_charge) {
      total+=time_q.total_seconds();
      time_q.zero_total();
    }
    if (_rot) {
      total+=time_quat.total_seconds();
      time_quat.zero_total();
    }
    if (_vel) {
      total+=time_vel.total_seconds();
      time_vel.zero_total();
    }

    return total+_time_transfer/1000.0;
  }

  /// Return the total time for data cast/pack
  /** Zeros the time so that atom times are only included once **/
  inline double cast_time()
    { double t=_time_cast; _time_cast=0.0; return t; }

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
    view.view_offset(0,buffer,m_size*m_size);
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
    view.view_offset(0,buffer,m_size*m_size);
    ucl_copy(dev_v,view,false);
  }

  /// Pack LAMMPS atom type constants into 2 vectors and copy to device
  template <class dev_typ, class t1, class t2>
  inline void type_pack2(const int n, UCL_D_Vec<dev_typ> &dev_v,
                         UCL_H_Vec<numtyp> &buffer, t1 ***one, t2 ***two) {
    int ii=0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        for (int k=0; k<n; k++) {
          buffer[ii*2]=static_cast<numtyp>(one[i][j][k]);
          buffer[ii*2+1]=static_cast<numtyp>(two[i][j][k]);
          ii++;
        }
      }
    }
    UCL_H_Vec<dev_typ> view;
    view.view_offset(0,buffer,n*n*n);
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
    view.view_offset(0,buffer,m_size*m_size);
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
    view.view_offset(0,buffer,m_size*m_size);
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
    view.view_offset(0,buffer,n);
    ucl_copy(dev_v,view,false);
  }

  // -------------------------COPY TO GPU ----------------------------------

  /// Signal that we need to transfer atom data for next timestep
  inline void data_unavail()
    { _x_avail=false; _q_avail=false; _quat_avail=false; _v_avail=false; _resized=false; }

  typedef struct { double x,y,z; } vec3d;
  typedef struct { numtyp x,y,z,w; } vec4d_t;

  /// Cast positions and types to write buffer
  inline void cast_x_data(double **host_ptr, const int *host_type) {
    if (_x_avail==false) {
      double t=MPI_Wtime();
      #ifdef GPU_CAST
      memcpy(host_x_cast.begin(),host_ptr[0],_nall*3*sizeof(double));
      memcpy(host_type_cast.begin(),host_type,_nall*sizeof(int));
      #else
      vec3d *host_p=reinterpret_cast<vec3d*>(&(host_ptr[0][0]));
      vec4d_t *xp=reinterpret_cast<vec4d_t*>(&(x[0]));
      #if (LAL_USE_OMP == 1)
      #pragma omp parallel for schedule(static)
      #endif
      for (int i=0; i<_nall; i++) {
        xp[i].x=host_p[i].x;
        xp[i].y=host_p[i].y;
        xp[i].z=host_p[i].z;
        xp[i].w=host_type[i];
      }
      #endif
      _time_cast+=MPI_Wtime()-t;
    }
  }

  /// Copy positions and types to device asynchronously
  /** Copies nall() elements **/
  inline void add_x_data(double **host_ptr, int *host_type) {
    time_pos.start();
    if (_x_avail==false) {
      #ifdef GPU_CAST
      x_cast.update_device(_nall*3,true);
      type_cast.update_device(_nall,true);
      int block_size=64;
      int GX=static_cast<int>(ceil(static_cast<double>(_nall)/block_size));
      k_cast_x.set_size(GX,block_size);
      k_cast_x.run(&x, &x_cast, &type_cast, &_nall);
      #else
      x.update_device(_nall*4,true);
      #endif
      _x_avail=true;
    }
    time_pos.stop();
  }

  /// Calls cast_x_data and add_x_data and times the routines
  inline void cast_copy_x(double **host_ptr, int *host_type) {
    cast_x_data(host_ptr,host_type);
    add_x_data(host_ptr,host_type);
  }

  // Cast charges to write buffer
  template<class cpytyp>
  inline void cast_q_data(cpytyp *host_ptr) {
    if (_q_avail==false) {
      double t=MPI_Wtime();
      // If double precision, still memcpy for async transfers
      if (_host_view) {
        q.host.view((numtyp*)host_ptr,_nall,*dev);
        q.device.view(q.host);
      } else if (sizeof(numtyp)==sizeof(double))
        memcpy(q.host.begin(),host_ptr,_nall*sizeof(numtyp));
      else
        #if (LAL_USE_OMP == 1) && (LAL_USE_OMP_SIMD == 1)
        #pragma omp parallel for simd schedule(static)
        #elif (LAL_USE_OMP_SIMD == 1)
        #pragma omp simd
        #endif
        for (int i=0; i<_nall; i++) q[i]=host_ptr[i];
      _time_cast+=MPI_Wtime()-t;
    }
  }

  // Copy charges to device asynchronously
  inline void add_q_data() {
    time_q.start();
    if (_q_avail==false) {
      q.update_device(_nall,true);
      _q_avail=true;
    }
    time_q.stop();
  }

  // Cast quaternions to write buffer
  template<class cpytyp>
  inline void cast_quat_data(cpytyp *host_ptr) {
    if (_quat_avail==false) {
      double t=MPI_Wtime();
      if (_host_view) {
        quat.host.view((numtyp*)host_ptr,_nall*4,*dev);
        quat.device.view(quat.host);
      } else if (sizeof(numtyp)==sizeof(double))
        memcpy(quat.host.begin(),host_ptr,_nall*4*sizeof(numtyp));
      else
        #if (LAL_USE_OMP == 1) && (LAL_USE_OMP_SIMD == 1)
        #pragma omp parallel for simd schedule(static)
        #elif (LAL_USE_OMP_SIMD == 1)
        #pragma omp simd
        #endif
        for (int i=0; i<_nall*4; i++) quat[i]=host_ptr[i];
      _time_cast+=MPI_Wtime()-t;
    }
  }

  // Copy quaternions to device
  /** Copies nall()*4 elements **/
  inline void add_quat_data() {
    time_quat.start();
    if (_quat_avail==false) {
      quat.update_device(_nall*4,true);
      _quat_avail=true;
    }
    time_quat.stop();
  }

  /// Cast velocities and tags to write buffer
  inline void cast_v_data(double **host_ptr, const tagint *host_tag) {
    if (_v_avail==false) {
      double t=MPI_Wtime();
      #ifdef GPU_CAST
      memcpy(host_v_cast.begin(),host_ptr[0],_nall*3*sizeof(double));
      memcpy(host_tag_cast.begin(),host_tag,_nall*sizeof(int));
      #else
      vec3d *host_p=reinterpret_cast<vec3d*>(&(host_ptr[0][0]));
      vec4d_t *vp=reinterpret_cast<vec4d_t*>(&(v[0]));
      #if (LAL_USE_OMP == 1)
      #pragma omp parallel for schedule(static)
      #endif
      for (int i=0; i<_nall; i++) {
        vp[i].x=host_p[i].x;
        vp[i].y=host_p[i].y;
        vp[i].z=host_p[i].z;
        vp[i].w=host_tag[i];
      }
      #endif
      _time_cast+=MPI_Wtime()-t;
    }
  }

  /// Copy velocities and tags to device asynchronously
  /** Copies nall() elements **/
  inline void add_v_data(double **host_ptr, tagint *host_tag) {
    time_vel.start();
    if (_v_avail==false) {
      #ifdef GPU_CAST
      v_cast.update_device(_nall*3,true);
      tag_cast.update_device(_nall,true);
      int block_size=64;
      int GX=static_cast<int>(ceil(static_cast<double>(_nall)/block_size));
      k_cast_x.set_size(GX,block_size);
      k_cast_x.run(&v, &v_cast, &tag_cast, &_nall);
      #else
      v.update_device(_nall*4,true);
      #endif
      _v_avail=true;
    }
    time_vel.stop();
  }

  /// Calls cast_v_data and add_v_data and times the routines
  inline void cast_copy_v(double **host_ptr, tagint *host_tag) {
    cast_v_data(host_ptr,host_tag);
    add_v_data(host_ptr,host_tag);
  }

  /// Add in casting time from additional data (seconds)
  inline void add_cast_time(double t) { _time_cast+=t; }

  /// Add in transfer time from additional data (ms)
  inline void add_transfer_time(double t) { _time_transfer+=t; }

  /// Return number of bytes used on device
  inline double max_gpu_bytes()
    { double m=_max_gpu_bytes; _max_gpu_bytes=0.0; return m; }

  /// Returns true if the device is addressing memory on the host
  inline bool host_view() { return _host_view; }

  // ------------------------------ DATA ----------------------------------

  /// Atom coordinates and types ([0] is x, [1] is y, [2] is z, [3] is type
  UCL_Vector<numtyp,numtyp> x;
  /// Charges
  UCL_Vector<numtyp,numtyp> q;
  /// Quaterions
  UCL_Vector<numtyp,numtyp> quat;
  /// Velocities
  UCL_Vector<numtyp,numtyp> v;

  #ifdef GPU_CAST
  UCL_Vector<double,double> x_cast;
  UCL_Vector<int,int> type_cast;
  #endif

  /// Cell list identifiers for device nbor builds
  UCL_D_Vec<unsigned> dev_cell_id;
  /// Cell list identifiers for device nbor builds
  UCL_D_Vec<int> dev_particle_id;

  /// Atom tag information for device nbor builds
  UCL_D_Vec<tagint> dev_tag;

  /// Cell list identifiers for hybrid nbor builds
  UCL_H_Vec<int> host_cell_id;
  /// Cell list identifiers for hybrid nbor builds
  UCL_H_Vec<int> host_particle_id;

  /// Device timers
  UCL_Timer time_pos, time_q, time_quat, time_vel;

  /// Geryon device
  UCL_Device *dev;

 private:
  #ifdef GPU_CAST
  UCL_Program *atom_program;
  UCL_Kernel k_cast_x;
  void compile_kernels(UCL_Device &dev);
  #endif

  bool _compiled;

  // True if data has been copied to device already
  bool _x_avail, _q_avail, _quat_avail, _v_avail, _resized;

  bool alloc(const int nall);

  bool _allocated, _rot, _charge, _bonds, _vel, _other;
  int _max_atoms, _nall, _gpu_nbor;
  bool _host_view;
  double _time_cast, _time_transfer;

  double _max_gpu_bytes;

  #ifdef USE_CUDPP
  CUDPPConfiguration sort_config;
  CUDPPHandle sort_plan;
  #endif

  #ifdef USE_HIP_DEVICE_SORT
  unsigned* sort_out_keys = nullptr;
  int* sort_out_values = nullptr;
  void* sort_temp_storage = nullptr;
  size_t sort_temp_storage_size = 0;
  size_t sort_out_size = 0;
  #endif
};

}

#endif

