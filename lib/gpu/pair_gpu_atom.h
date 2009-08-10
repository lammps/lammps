/***************************************************************************
                               pair_gpu_atom.h
                             -------------------
                               W. Michael Brown

  Memory routines for moving atom and force data between host and gpu

 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Tue Aug 4 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#ifndef PAIR_GPU_ATOM_H
#define PAIR_GPU_ATOM_H

// PRECISION - Precision for rsq, energy, force, and torque calculation
// ACC_PRECISION - Precision for accumulation of energies, forces, and torques
#ifdef _SINGLE_DOUBLE
#define PRECISION float
#define ACC_PRECISION double
#endif

#ifdef _DOUBLE_DOUBLE
#define PRECISION double
#define ACC_PRECISION double
#endif

#ifndef PRECISION
#define PRECISION float
#define ACC_PRECISION double
#endif

#define MAX_ATOMS 65536
#include "nvc_timer.h"
#include "pair_gpu_texture.h"

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
  void init(const int max_atoms);
  
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
  
  /// Copy num_rows x nall() write buffer to x in GPU
  /** num_rows<=atom_fields() **/
  inline void copy_atom_data(const int num_rows, cudaStream_t &stream) 
    { dev_x.copy_2Dfrom_host(host_write.begin(),num_rows,nall(),stream); }

    
  // -------------------------COPY FROM GPU -------------------------------

  /// Copy answers from GPU into read buffer
  void copy_answers(const bool eflag, const bool vflag, cudaStream_t &s);
  
  /// Copy energy and virial data into LAMMPS memory
  double energy_virial(const int *ilist, const bool eflag_atom,
                       const bool vflag_atom, double *eatom, double **vatom,
                       double *virial);

  /// Add forces from the GPU into a LAMMPS pointer
  void add_forces(const int *ilist, double **f);

  /// Add torques from the GPU into a LAMMPS pointer
  void add_torques(const int *ilist, double **tor, const int n);
    
  // ------------------------------ DATA ----------------------------------

  // atom_fields() x n (example: rows 1-3 position, 4 is type)
  NVC_ConstMatT dev_x;
  // ans_fields() x n Storage for Forces, etc.
  // example: if (eflag and vflag) 1 is energy, 2-7 is virial, 8-10 is force 
  // example: if (!eflag) 1-3 is force
  NVC_Mat<acctyp> ans;                               

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
