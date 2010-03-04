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
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#ifndef PAIR_GPU_ATOM_H
#define PAIR_GPU_ATOM_H

// PRECISION - Precision for rsq, energy, force, and torque calculation
// ACC_PRECISION - Precision for accumulation of energies, forces, and torques
#ifdef _SINGLE_DOUBLE
#define PRECISION float
#define ACC_PRECISION double
#define MAX_ATOMS 65536
#define vec4 float4
#endif

#ifdef _DOUBLE_DOUBLE
#define PRECISION double
#define ACC_PRECISION double
#define MAX_ATOMS 32768
struct vec4 { double x; double y; double z; double w; };
#endif

#ifndef PRECISION
#define PRECISION float
#define ACC_PRECISION float
#define MAX_ATOMS 65536
#define vec4 float4
#endif

#include "nvc_timer.h"
#include "nvc_memory.h"

template <class numtyp, class acctyp>
class PairGPUAtom {
 public:
  PairGPUAtom() : _atom_fields(4), _ans_fields(10), allocated(false) {}
  ~PairGPUAtom() { clear(); }

  // Accessors
  inline int atom_fields() const { return _atom_fields; }
  inline int ans_fields() const { return _ans_fields; }
  inline int max_atoms() const { return _max_atoms; }
  inline int nall() const { return _nall; }
  inline int inum() const { return _inum; }

  /// Set number of atoms for future copy operations
  inline void nall(const int n) { _nall=n; }
  /// Set number of inum for future copy operations
  inline void inum(const int n) { _inum=n; }
  /// Set the number of atom fields (x, y, z, type, etc)
  inline void atom_fields(const int n) { _atom_fields=n; }
  /// Set the number of answer fields (energy, virial, force, etc.)
  inline void ans_fields(const int n) { _ans_fields=n; }
  
  /// Memory usage per atom in this class
  /** \note atom_fields and ans_fields should be set for correct answer **/
  int bytes_per_atom() const; 

  /// Must be called once to allocate host and device memory
  /** \note atom_fields and ans_fields should be set first if not default **/
  bool init(const int max_atoms);
  void resize(const int max_atoms, bool &success);
  
  /// Free all memory on host and device
  void clear();
 
  /// Return the total amount of host memory used by class
  double host_memory_usage(const int max_atoms) const;

  
  // -------------------------COPY TO GPU ----------------------------------

  /// Reset the write buffer pointer (Start copying new atom data)
  inline void reset_write_buffer() { _write_loc=host_write.begin(); }
  
  /// Add a row to write buffer with unit stride
  /** Copies nall() elements **/
  template<class cpytyp>
  inline void add_atom_data(const cpytyp *host_ptr)
    { for (int i=0; i<_nall; i++) { *_write_loc=host_ptr[i]; _write_loc++; } }
  
  /// Add a row to write buffer with non-unit stride
  /** Copies nall() elements **/
  template<class cpytyp>
  inline void add_atom_data(const cpytyp *hostptr, const int stride) {
    int t=_nall*stride; 
    for (int i=0; i<t; i+=stride) { *_write_loc=hostptr[i]; _write_loc++; }
  }
  
  /// Add positions to write buffer
  /** Copies nall() elements **/
  inline void add_x_data(double **host_ptr, const int *host_type) {
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
  }      

  /// Add quaternions to write buffer
  /** Copies nall() elements **/
  template<class cpytyp>
  inline void add_q_data(const cpytyp *host_ptr) {
    const int end=_nall*4;
    for (int i=0; i<end; i++) { *_write_loc=host_ptr[i]; _write_loc++; } 
  } 

  /// Copy num_rows positions+type to x in GPU
  /** num_rows<=atom_fields() **/
  inline void copy_x_data(cudaStream_t &stream) 
    { dev_x.copy_from_host(host_write.begin(),_nall*4,stream); }
  inline void copy_q_data(cudaStream_t &stream) 
    { dev_q.copy_from_host(host_write.begin()+_nall*4,_nall*4,stream); }
    
  // -------------------------COPY FROM GPU -------------------------------

  /// Copy answers from GPU into read buffer
  void copy_answers(const bool eflag, const bool vflag, cudaStream_t &s);
  
  /// Copy energy and virial data into LAMMPS memory
  double energy_virial(const int *ilist, const bool eflag_atom,
                       const bool vflag_atom, double *eatom, double **vatom,
                       double *virial, double **f, double **tor, const int);
                       
  /// Add forces and torques from the GPU into a LAMMPS pointer
  void copy_asphere(const int *ilist, double **f, double **tor, const int n);
  // ------------------------------ DATA ----------------------------------

  // atom coordinates
  NVC_Vec<numtyp> dev_x;
  // quaterions
  NVC_Vec<numtyp> dev_q;
  // ans_fields()
  // example: if (eflag and vflag) 1 is energy, 2-7 is virial
  NVC_Vec<acctyp> ans;                               

  // Buffer for moving floating point data to GPU
  NVC_HostT host_write;
  // Buffer for moving floating point data to CPU
  NVC_Host<acctyp> host_read;
  
  // Timing Stuff
  NVCTimer time_atom, time_answer;
  
 private:
  bool allocated, _eflag, _vflag;
  int _atom_fields, _ans_fields;
  int _max_atoms, _nall, _inum;
  numtyp * _write_loc;
  acctyp * _read_loc;
};

#endif
