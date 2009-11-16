/***************************************************************************
                                  lj_gpu.cu
                             -------------------
                               W. Michael Brown

  Lennard-Jones potential GPU calcultation

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

#include <iostream>
#include <cassert>
#include "nvc_macros.h"
#include "nvc_timer.h"
#include "nvc_device.h"
#include "pair_gpu_texture.h"
#include "lj_gpu_memory.cu"
#include "lj_gpu_kernel.h"

using namespace std;

static LJ_GPU_Memory<PRECISION,ACC_PRECISION> LJMF;
#define LJMT LJ_GPU_Memory<numtyp,acctyp>

// ---------------------------------------------------------------------------
// Convert something to a string
// ---------------------------------------------------------------------------
#include <sstream>
template <class t>
inline string lj_gpu_toa(const t& in) {
  ostringstream o;
  o.precision(2);
  o << in;
  return o.str();
}

// ---------------------------------------------------------------------------
// Return string with GPU info
// ---------------------------------------------------------------------------
string lj_gpu_name(const int id, const int max_nbors) {
  string name=LJMF.gpu.name(id)+", "+
              lj_gpu_toa(LJMF.gpu.cores(id))+" cores, "+
              lj_gpu_toa(LJMF.gpu.gigabytes(id))+" GB, "+
              lj_gpu_toa(LJMF.gpu.clock_rate(id))+" GHZ, "+
              lj_gpu_toa(LJMF.get_max_atoms(LJMF.gpu.bytes(id),
                                            max_nbors))+" Atoms";
  return name;
}

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
bool lj_gpu_init(int &ij_size, const int ntypes, double **cutsq,double **sigma, 
                 double **epsilon, double **host_lj1, double **host_lj2, 
                 double **host_lj3, double **host_lj4, double **offset, 
                 double *special_lj, const int max_nbors, const int gpu_id) {
  LJMF.gpu.init();                  
  if (LJMF.gpu.num_devices()==0)
    return false;                   

  ij_size=IJ_SIZE;
  return LJMF.init(ij_size, ntypes, cutsq, sigma, epsilon, host_lj1, host_lj2, 
                   host_lj3, host_lj4, offset, special_lj, max_nbors, gpu_id);
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
void lj_gpu_clear() {
  LJMF.clear();
}

// ---------------------------------------------------------------------------
// copy atom positions and optionally types to device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline void _lj_gpu_atom(PairGPUAtom<numtyp,acctyp> &atom, double **host_x,
                         const int *host_type, const bool rebuild,
                         cudaStream_t &stream) {
  atom.time_atom.start();
  atom.reset_write_buffer();
  
  // First row of dev_x is x position, second is y, third is z
  atom.add_atom_data(host_x[0],3);
  atom.add_atom_data(host_x[0]+1,3);
  atom.add_atom_data(host_x[0]+2,3);

  int csize=3;

  // If a rebuild occured, copy type data
  if (rebuild) {
    atom.add_atom_data(host_type);
    csize++;
  }
  
  atom.copy_atom_data(csize,stream);
  atom.time_atom.stop();
}

void lj_gpu_atom(double **host_x, const int *host_type, const bool rebuild) {
  _lj_gpu_atom(LJMF.atom, host_x, host_type, rebuild, LJMF.pair_stream);
}

// ---------------------------------------------------------------------------
// Signal that we need to transfer a new neighbor list
// ---------------------------------------------------------------------------
template <class LJMTyp>
bool _lj_gpu_reset_nbors(LJMTyp &ljm, const int nall, const int inum, 
                         int *ilist, const int *numj) {
  if (nall>ljm.max_atoms)
    return false;
  
  ljm.nbor.time_nbor.start();

  ljm.atom.nall(nall);
  ljm.atom.inum(inum);
  ljm.nbor.reset(inum,ilist,numj,ljm.pair_stream);

  ljm.nbor.time_nbor.stop();
  return true;
}

bool lj_gpu_reset_nbors(const int nall, const int inum, int *ilist, 
                        const int *numj) {
  return _lj_gpu_reset_nbors(LJMF,nall,inum,ilist,numj);
}

// ---------------------------------------------------------------------------
// Copy a set of ij_size ij interactions to device and compute energies,
// forces, and torques for those interactions
// ---------------------------------------------------------------------------
template <class LJMTyp>
void _lj_gpu_nbors(LJMTyp &ljm, const int *ij, const int num_ij) {
  ljm.nbor.time_nbor.add_to_total();
  
  // CUDA_SAFE_CALL(cudaStreamSynchronize(ljm.pair_stream)); // Not if timed

  memcpy(ljm.nbor.host_ij.begin(),ij,num_ij*sizeof(int));
  ljm.nbor.time_nbor.start();
  ljm.nbor.add(num_ij,ljm.pair_stream);
  ljm.nbor.time_nbor.stop();
}

void lj_gpu_nbors(const int *ij, const int num_ij) {
  _lj_gpu_nbors(LJMF,ij,num_ij);
}

// ---------------------------------------------------------------------------
// Calculate energies and forces for all ij interactions
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void _lj_gpu(LJMT &ljm, const bool eflag, const bool vflag, const bool rebuild){
  // Compute the block size and grid size to keep all cores busy
  const int BX=BLOCK_1D;

  int GX=static_cast<int>(ceil(static_cast<double>(ljm.atom.inum())/BX));

  ljm.time_pair.start();
  if (ljm.shared_types)
    kernel_lj_fast<numtyp,acctyp><<<GX,BX,0,ljm.pair_stream>>>
           (ljm.special_lj.begin(), ljm.nbor.dev_nbor.begin(), 
            ljm.nbor.ij.begin(), ljm.nbor.dev_nbor.row_size(), 
            ljm.atom.ans.begin(), ljm.atom.ans.row_size(), eflag,
            vflag, ljm.atom.inum(), ljm.atom.nall());
  else
    kernel_lj<numtyp,acctyp><<<GX,BX,0,ljm.pair_stream>>>
           (ljm.special_lj.begin(), ljm.nbor.dev_nbor.begin(), 
            ljm.nbor.ij.begin(), ljm.nbor.dev_nbor.row_size(), 
            ljm.atom.ans.begin(), ljm.atom.ans.row_size(), eflag, 
            vflag, ljm.atom.inum(), ljm.atom.nall());
  ljm.time_pair.stop();
}

void lj_gpu(const bool eflag, const bool vflag, const bool rebuild) {
  _lj_gpu<PRECISION,ACC_PRECISION>(LJMF,eflag,vflag,rebuild);
}

// ---------------------------------------------------------------------------
// Get energies and forces to host
// ---------------------------------------------------------------------------
template<class numtyp, class acctyp>
double _lj_gpu_forces(LJMT &ljm, double **f, const int *ilist,
                      const bool eflag, const bool vflag, const bool eflag_atom,
                      const bool vflag_atom, double *eatom, double **vatom,
                      double *virial) {
  double evdw;
  
  ljm.atom.time_answer.start();
  ljm.atom.copy_answers(eflag,vflag,ljm.pair_stream);

  ljm.atom.time_atom.add_to_total();
  ljm.nbor.time_nbor.add_to_total();
  ljm.time_pair.add_to_total();
  CUDA_SAFE_CALL(cudaStreamSynchronize(ljm.pair_stream));

  evdw=ljm.atom.energy_virial(ilist,eflag_atom,vflag_atom,eatom,vatom,virial);
  ljm.atom.add_forces(ilist,f);
  ljm.atom.time_answer.stop();
  ljm.atom.time_answer.add_to_total();
  return evdw;
}

double lj_gpu_forces(double **f, const int *ilist, const bool eflag, 
                     const bool vflag, const bool eflag_atom,
                     const bool vflag_atom, double *eatom, double **vatom,
                     double *virial) {
  return _lj_gpu_forces<PRECISION,ACC_PRECISION> 
    (LJMF,f,ilist,eflag,vflag,eflag_atom,vflag_atom,eatom,vatom,virial);
}

void lj_gpu_time() {
  cout.precision(4);
  cout << "Atom copy:     " << LJMF.atom.time_atom.total_seconds() << " s.\n";
  cout << "Neighbor copy: " << LJMF.nbor.time_nbor.total_seconds() << " s.\n";
  cout << "LJ calc:       " << LJMF.time_pair.total_seconds() << " s.\n";
  cout << "Answer copy:   " << LJMF.atom.time_answer.total_seconds() << " s.\n";
}

int lj_gpu_num_devices() {
  return LJMF.gpu.num_devices();
}

double lj_gpu_bytes() {
  return LJMF.host_memory_usage();
}
