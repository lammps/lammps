/***************************************************************************
                                  gb_gpu.cu
                             -------------------
                               W. Michael Brown

  Gay-Berne anisotropic potential GPU calcultation

   *** Force decomposition by Atom Version ***

 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Tue Jun 23 2009
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
#include "gb_gpu_memory.cu"
#include "gb_gpu_kernel.h"

using namespace std;

static GB_GPU_Memory<PRECISION,ACC_PRECISION> GBMF[MAX_GPU_THREADS];
#define GBMT GB_GPU_Memory<numtyp,acctyp>

// ---------------------------------------------------------------------------
// Pack neighbors from dev_ij array into dev_nbor matrix for coalesced access
// -- Only pack neighbors matching the specified inclusive range of forms
// -- Only pack neighbors within cutoff
// ---------------------------------------------------------------------------
template<class numtyp>
__global__ void kernel_pack_nbor(int *dev_nbor, const int nbor_pitch, 
                                 const int start, const int inum, 
                                 const int *dev_ij, const int form_low, 
                                 const int form_high, const int nall) {
                                
  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x+INT_MUL(blockIdx.x,blockDim.x)+start;

  if (ii<inum) {
    int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
    nbor+=nbor_pitch;
    int *nbor_newj=nbor;
    nbor+=nbor_pitch;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=_x_<numtyp>(i,7);

    int newj=0;  
    for ( ; list<list_end; list++) {
      int j=*list;
      if (j>=nall)
        j%=nall;
      int jtype=_x_<numtyp>(j,7);
      
      if (_form_(itype,jtype)>=form_low && _form_(itype,jtype)<=form_high) {
        // Compute r12;
        numtyp rsq=_x_<numtyp>(j,0)-ix;
        rsq*=rsq;
        numtyp t=_x_<numtyp>(j,1)-iy;
        rsq+=t*t;
        t=_x_<numtyp>(j,2)-iz;
        rsq+=t*t;

        if (rsq< _cutsq_<numtyp>(itype,jtype)) {
          *nbor=j;
          nbor+=nbor_pitch;
          newj++;
        }
      }
    }
    *nbor_newj=newj;
  }
}

// ---------------------------------------------------------------------------
// Pack neighbors from dev_ij array into dev_nbor matrix for coalesced access
// -- Only pack neighbors matching the specified inclusive range of forms
// -- Only pack neighbors within cutoff
// -- Fast version of routine that uses shared memory for LJ constants
// ---------------------------------------------------------------------------
template<class numtyp>
__global__ void kernel_pack_nbor_fast(int *dev_nbor, const int nbor_pitch, 
                                      const int start, const int inum, 
                                      const int *dev_ij, const int form_low, 
                                      const int form_high, const int nall) {
                                
  int ii=threadIdx.x;
  __shared__ int form[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (ii<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    int itype=ii/MAX_SHARED_TYPES;
    int jtype=ii%MAX_SHARED_TYPES;
    cutsq[ii]=_cutsq_<numtyp>(itype,jtype);
    form[ii]=_form_(itype,jtype);
  }
  ii+=INT_MUL(blockIdx.x,blockDim.x)+start;
  __syncthreads();

  if (ii<inum) {
    int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
    nbor+=nbor_pitch;
    int *nbor_newj=nbor;
    nbor+=nbor_pitch;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=INT_MUL(MAX_SHARED_TYPES,_x_<numtyp>(i,7));

    int newj=0;  
    for ( ; list<list_end; list++) {
      int j=*list;
      if (j>=nall)
        j%=nall;
      int mtype=itype+_x_<numtyp>(j,7);
      
      if (form[mtype]>=form_low && form[mtype]<=form_high) {
        // Compute r12;
        numtyp rsq=_x_<numtyp>(j,0)-ix;
        rsq*=rsq;
        numtyp t=_x_<numtyp>(j,1)-iy;
        rsq+=t*t;
        t=_x_<numtyp>(j,2)-iz;
        rsq+=t*t;

        if (rsq<cutsq[mtype]) {
          *nbor=j;
          nbor+=nbor_pitch;
          newj++;
        }
      }
    }
    *nbor_newj=newj;
  }
}

template<class numtyp, class acctyp>
void pack_nbors(GBMT &gbm, const int GX, const int BX, const int start, 
                const int inum, const int form_low, const int form_high) {
  if (gbm.shared_types)
    kernel_pack_nbor_fast<numtyp><<<GX,BX,0,gbm.pair_stream>>>
          (gbm.nbor.dev_nbor.begin(), gbm.nbor.dev_nbor.row_size(), start, inum,
           gbm.nbor.ij.begin(),form_low,form_high,gbm.atom.nall());
  else
    kernel_pack_nbor<numtyp><<<GX,BX,0,gbm.pair_stream>>>
          (gbm.nbor.dev_nbor.begin(), gbm.nbor.dev_nbor.row_size(), start, inum,
           gbm.nbor.ij.begin(),form_low,form_high,gbm.atom.nall());
}

// ---------------------------------------------------------------------------
// Convert something to a string
// ---------------------------------------------------------------------------
#include <sstream>
template <class t>
inline string gb_gpu_toa(const t& in) {
  ostringstream o;
  o.precision(2);
  o << in;
  return o.str();
}

// ---------------------------------------------------------------------------
// Return string with GPU info
// ---------------------------------------------------------------------------
string gb_gpu_name(const int id, const int max_nbors) {
  string name=GBMF[0].gpu.name(id)+", "+
              gb_gpu_toa(GBMF[0].gpu.cores(id))+" cores, "+
              gb_gpu_toa(GBMF[0].gpu.gigabytes(id))+" GB, "+
              gb_gpu_toa(GBMF[0].gpu.clock_rate(id))+" GHZ, "+
              gb_gpu_toa(GBMF[0].get_max_atoms(GBMF[0].gpu.bytes(id),
                                               max_nbors))+" Atoms";
  return name;
}

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int * gb_gpu_init(int &ij_size, const int ntypes, const double gamma,
                  const double upsilon, const double mu, double **shape,
                  double **well, double **cutsq, double **sigma, 
                  double **epsilon, double *host_lshape, int **form,
                  double **host_lj1, double **host_lj2, double **host_lj3, 
                  double **host_lj4, double **offset, double *special_lj,
                  const int max_nbors, const int thread, const int gpu_id) {
  assert(thread<MAX_GPU_THREADS);
  
  if (GBMF[thread].gpu.num_devices()==0)
    return 0;                   

  ij_size=IJ_SIZE;
  return GBMF[thread].init(ij_size, ntypes, gamma, upsilon, mu, shape,
                           well, cutsq, sigma, epsilon, host_lshape, form,
                           host_lj1, host_lj2, host_lj3, host_lj4, offset,
                           special_lj, max_nbors, false, gpu_id);
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
void gb_gpu_clear(const int thread) {
  GBMF[thread].clear();
}

// ---------------------------------------------------------------------------
// copy atom positions, quaternions, and optionally types to device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline void _gb_gpu_atom(PairGPUAtom<numtyp,acctyp> &atom, double **host_x, 
                         double **host_quat, const int *host_type, 
                         const bool rebuild, cudaStream_t &stream) {
  atom.time_atom.start();
  atom.reset_write_buffer();
 
  // Rows 1-3 of dev_x are position; rows 4-7 are quaternion
  atom.add_atom_data(host_x[0],3);
  atom.add_atom_data(host_x[0]+1,3);
  atom.add_atom_data(host_x[0]+2,3);
  atom.add_atom_data(host_quat[0],4);
  atom.add_atom_data(host_quat[0]+1,4);
  atom.add_atom_data(host_quat[0]+2,4);
  atom.add_atom_data(host_quat[0]+3,4);
   
  int csize=7;
  
  // If a rebuild occured, copy type data
  if (rebuild) {
    atom.add_atom_data(host_type);
    csize++;
  }
  
  atom.copy_atom_data(csize,stream);
  atom.time_atom.stop();
}

void gb_gpu_atom(double **host_x, double **host_quat, 
                 const int *host_type, const bool rebuild, const int thread) {
  _gb_gpu_atom(GBMF[thread].atom, host_x, host_quat, host_type, rebuild,
               GBMF[thread].pair_stream);
}

// ---------------------------------------------------------------------------
// Signal that we need to transfer a new neighbor list
// ---------------------------------------------------------------------------
template <class gbmtyp>
int * _gb_gpu_reset_nbors(gbmtyp &gbm, const int nall, const int nlocal, 
                          const int inum, int *ilist, const int *numj,
                          const int *type, bool &success) {
  if (nall>gbm.max_atoms) {
    success=false;
    return 0;
  }
  success=true;
    
  gbm.nbor.time_nbor.start();

  gbm.atom.nall(nall);
  gbm.atom.inum(inum);

  if (gbm.multiple_forms) {
    int p=0, acc=0;
    for (int i=0; i<inum; i++) {
      int itype=type[ilist[i]];
      if (gbm.host_form[itype][itype]==ELLIPSE_ELLIPSE) {
        gbm.host_olist[p]=ilist[i];
        gbm.nbor.host_ij[p]=numj[ilist[i]];
        gbm.nbor.host_ij[p+inum]=acc;
        acc+=numj[ilist[i]];
        p++;
      }
    }
    gbm.last_ellipse=p;
    for (int i=0; i<inum; i++) {
      int itype=type[ilist[i]];
      if (gbm.host_form[itype][itype]!=ELLIPSE_ELLIPSE) {
        gbm.host_olist[p]=ilist[i];
        gbm.nbor.host_ij[p]=numj[ilist[i]];
        gbm.nbor.host_ij[p+inum]=acc;
        acc+=numj[ilist[i]];
        p++;
      }
    }
    gbm.nbor.ij_total=0;
    gbm.nbor.dev_nbor.copy_from_host(gbm.host_olist.begin(),inum);
    gbm.nbor.host_ij.copy_to_2Ddevice(gbm.nbor.dev_nbor.begin()+
                                      gbm.nbor.dev_nbor.row_size(),
                                      gbm.nbor.dev_nbor.row_size(),2,inum,
                                      gbm.pair_stream);
  } else {
    gbm.nbor.reset(inum,ilist,numj,gbm.pair_stream);
    gbm.last_ellipse=inum;
  }

  gbm.nbor.time_nbor.stop();
  
  if (gbm.multiple_forms)
    return gbm.host_olist.begin();
  return ilist;
}

int * gb_gpu_reset_nbors(const int nall, const int nlocal, const int inum, 
                         int *ilist, const int *numj, const int *type,
                         const int thread, bool &success) {
  return _gb_gpu_reset_nbors(GBMF[thread],nall,nlocal,inum,ilist,numj,type,
                             success);
}

// ---------------------------------------------------------------------------
// Copy a set of ij_size ij interactions to device and compute energies,
// forces, and torques for those interactions
// ---------------------------------------------------------------------------
template <class gbmtyp>
void _gb_gpu_nbors(gbmtyp &gbm, const int num_ij, const bool eflag) {
  gbm.nbor.time_nbor.add_to_total();
  // CUDA_SAFE_CALL(cudaStreamSynchronize(gbm.pair_stream)); // Not if timed
  
  gbm.nbor.time_nbor.start();
  gbm.nbor.add(num_ij,gbm.pair_stream);
  gbm.nbor.time_nbor.stop();
}

void gb_gpu_nbors(const int num_ij, const bool eflag, const int thread) {
  _gb_gpu_nbors(GBMF[thread],num_ij,eflag);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques for all ij interactions
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void _gb_gpu_gayberne(GBMT &gbm, const bool eflag, const bool vflag, 
                      const bool rebuild) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=BLOCK_1D;

  int GX=static_cast<int>(ceil(static_cast<double>(gbm.atom.inum())/BX));

  if (gbm.multiple_forms) {
    gbm.time_kernel.start();
    if (gbm.last_ellipse>0) {
      // ------------ ELLIPSE_ELLIPSE and ELLIPSE_SPHERE ---------------
      GX=static_cast<int>(ceil(static_cast<double>(gbm.last_ellipse)/
                               static_cast<double>(BX)));
      pack_nbors(gbm,GX,BX, 0, gbm.last_ellipse,SPHERE_ELLIPSE,ELLIPSE_ELLIPSE);
      gbm.time_kernel.stop();
  
      gbm.time_gayberne.start();                                 
      kernel_gayberne<numtyp,acctyp><<<GX,BX,0,gbm.pair_stream>>>
           (gbm.gamma_upsilon_mu.begin(), gbm.special_lj.begin(), 
            gbm.nbor.dev_nbor.begin(), gbm.nbor.dev_nbor.row_size(),
            gbm.atom.ans.begin(), gbm.atom.ans.row_size(),gbm.dev_error.begin(),
            eflag, vflag, gbm.last_ellipse, gbm.atom.nall());
      gbm.time_gayberne.stop();

      if (gbm.last_ellipse==gbm.atom.inum()) {
        gbm.time_kernel2.start();
        gbm.time_kernel2.stop();
        gbm.time_gayberne2.start();
        gbm.time_gayberne2.stop();
        gbm.time_pair.start();
        gbm.time_pair.stop();
        return;
      }
                    
      // ------------ SPHERE_ELLIPSE ---------------
        
      gbm.time_kernel2.start();
      GX=static_cast<int>(ceil(static_cast<double>(gbm.atom.inum()-
                               gbm.last_ellipse)/BX));
      pack_nbors(gbm,GX,BX,gbm.last_ellipse,gbm.atom.inum(),ELLIPSE_SPHERE,
                 ELLIPSE_SPHERE);
      gbm.time_kernel2.stop();

      gbm.time_gayberne2.start();
      kernel_sphere_gb<numtyp,acctyp><<<GX,BX,0,gbm.pair_stream>>>
           (gbm.gamma_upsilon_mu.begin(), gbm.special_lj.begin(), 
            gbm.nbor.dev_nbor.begin(), gbm.nbor.dev_nbor.row_size(), 
            gbm.atom.ans.begin(), gbm.atom.ans.row_size(),gbm.dev_error.begin(),
            eflag, vflag, gbm.last_ellipse, gbm.atom.inum(), gbm.atom.nall());
      gbm.time_gayberne2.stop();
   } else {
      gbm.atom.ans.zero();
      gbm.time_kernel.stop();
      gbm.time_gayberne.start();                                 
      gbm.time_gayberne.stop();
      gbm.time_kernel2.start();
      gbm.time_kernel2.stop();
      gbm.time_gayberne2.start();
      gbm.time_gayberne2.stop();
    }
    
    // ------------         LJ      ---------------
    gbm.time_pair.start();
    if (gbm.last_ellipse<gbm.atom.inum()) {
      if (gbm.shared_types)
        kernel_lj_fast<numtyp,acctyp><<<GX,BX,0,gbm.pair_stream>>>
           (gbm.special_lj.begin(), gbm.nbor.dev_nbor.begin(), 
            gbm.nbor.ij.begin(), gbm.nbor.dev_nbor.row_size(), 
            gbm.atom.ans.begin(), gbm.atom.ans.row_size(), 
            gbm.dev_error.begin(), eflag, vflag, gbm.last_ellipse, 
            gbm.atom.inum(), gbm.atom.nall());
      else
        kernel_lj<numtyp,acctyp><<<GX,BX,0,gbm.pair_stream>>>
           (gbm.special_lj.begin(), gbm.nbor.dev_nbor.begin(), 
            gbm.nbor.ij.begin(), gbm.nbor.dev_nbor.row_size(), 
            gbm.atom.ans.begin(), gbm.atom.ans.row_size(),gbm.dev_error.begin(),
            eflag, vflag, gbm.last_ellipse, gbm.atom.inum(), gbm.atom.nall());
    }
    gbm.time_pair.stop();
  } else {
    gbm.time_kernel.start();
    pack_nbors(gbm, GX, BX, 0, gbm.atom.inum(),SPHERE_SPHERE,ELLIPSE_ELLIPSE);
    gbm.time_kernel.stop();
  
    gbm.time_gayberne.start();                                 
    kernel_gayberne<numtyp,acctyp><<<GX,BX,0,gbm.pair_stream>>>
         (gbm.gamma_upsilon_mu.begin(), gbm.special_lj.begin(), 
          gbm.nbor.dev_nbor.begin(), gbm.nbor.dev_nbor.row_size(), 
          gbm.atom.ans.begin(), gbm.atom.ans.row_size(), gbm.dev_error.begin(), 
          eflag, vflag, gbm.atom.inum(), gbm.atom.nall());
    gbm.time_gayberne.stop();
  }
}

void gb_gpu_gayberne(const bool eflag, const bool vflag, const bool rebuild, 
                     const int thread) {
  _gb_gpu_gayberne<PRECISION,ACC_PRECISION>(GBMF[thread],eflag,vflag,rebuild);
}

// ---------------------------------------------------------------------------
// Get energies, forces, and torques to host
// ---------------------------------------------------------------------------
template<class numtyp, class acctyp>
double _gb_gpu_forces(GBMT &gbm, double **f, double **tor, const int *ilist,
                      const bool eflag, const bool vflag, const bool eflag_atom,
                      const bool vflag_atom, double *eatom, double **vatom,
                      double *virial) {
  double evdw;

  gbm.atom.time_answer.start();
  gbm.atom.copy_answers(eflag,vflag,gbm.pair_stream);

  gbm.atom.time_atom.add_to_total();
  gbm.nbor.time_nbor.add_to_total();
  gbm.time_kernel.add_to_total();
  gbm.time_gayberne.add_to_total();
  if (gbm.multiple_forms) {
    gbm.time_kernel2.add_to_total();
    gbm.time_gayberne2.add_to_total();
    gbm.time_pair.add_to_total();
  }      
  // CUDA_SAFE_CALL(cudaStreamSynchronize(gbm.pair_stream)); // Not if timed
  
  evdw=gbm.atom.energy_virial(ilist,eflag_atom,vflag_atom,eatom,vatom,virial);
  gbm.atom.add_forces(ilist,f);
  gbm.atom.add_torques(ilist,tor,gbm.last_ellipse);
  gbm.atom.time_answer.stop();
  gbm.atom.time_answer.add_to_total();
  return evdw;
}

double gb_gpu_forces(double **f, double **tor, const int *ilist,
                     const bool eflag, const bool vflag, const bool eflag_atom,
                     const bool vflag_atom, double *eatom, double **vatom,
                     double *virial, const int thread) {
  return _gb_gpu_forces<PRECISION,ACC_PRECISION>
                       (GBMF[thread],f,tor,ilist,eflag,vflag,eflag_atom,
                        vflag_atom,eatom,vatom,virial);
}

void gb_gpu_time(const int i) {
  cout.precision(4);
  cout << "Atom copy:     " << GBMF[i].atom.time_atom.total_seconds() 
       << " s.\n"
       << "Neighbor copy: " << GBMF[i].nbor.time_nbor.total_seconds() 
       << " s.\n"
       << "Neighbor pack: " << GBMF[i].time_kernel.total_seconds()+
                               GBMF[i].time_kernel2.total_seconds() << " s.\n"
       << "Force calc:    " << GBMF[i].time_gayberne.total_seconds()+
                               GBMF[i].time_gayberne2.total_seconds()<< " s.\n";
  if (GBMF[i].multiple_forms)
    cout << "LJ calc:       " << GBMF[i].time_pair.total_seconds() << " s.\n";
  cout << "Answer copy:   " << GBMF[i].atom.time_answer.total_seconds() 
       << " s.\n";
}

int gb_gpu_num_devices() {
  return GBMF[0].gpu.num_devices();
}

double gb_gpu_bytes() {
  return GBMF[0].host_memory_usage();
}
