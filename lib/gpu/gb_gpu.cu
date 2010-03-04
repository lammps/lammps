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
__global__ void kernel_pack_nbor(const vec4 *x_, int *dev_nbor, const int nbor_pitch, 
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
  
    vec4 ix=x_[i];
    int itype=ix.w;

    int newj=0;  
    for ( ; list<list_end; list++) {
      int j=*list;
      if (j>=nall)
        j%=nall;
      vec4 jx=x_[j];
      int jtype=jx.w;
      
      if (_form_(itype,jtype)>=form_low && _form_(itype,jtype)<=form_high) {
        // Compute r12;
        numtyp rsq=jx.x-ix.x;
        rsq*=rsq;
        numtyp t=jx.y-ix.y;
        rsq+=t*t;
        t=jx.z-ix.z;
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
__global__ void kernel_pack_nbor_fast(const vec4 *x_, int *dev_nbor, const int nbor_pitch, 
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
  
    vec4 ix=x_[i];
    int itype=INT_MUL(MAX_SHARED_TYPES,ix.w);

    int newj=0;  
    for ( ; list<list_end; list++) {
      int j=*list;
      if (j>=nall)
        j%=nall;
      vec4 jx=x_[j];
      int jtype=jx.w;
      int mtype=itype+jtype;
      
      if (form[mtype]>=form_low && form[mtype]<=form_high) {
        // Compute r12;
        numtyp rsq=jx.x-ix.x;
        rsq*=rsq;
        numtyp t=jx.y-ix.y;
        rsq+=t*t;
        t=jx.z-ix.z;
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
  if (gbm.shared_types) {
    kernel_pack_nbor_fast<numtyp><<<GX,BX,0,gbm.pair_stream>>>
          ((vec4 *)gbm.atom.dev_x.begin(),gbm.nbor.dev_nbor.begin(), 
           gbm.atom.inum(), start, inum,
           gbm.nbor.ij.begin(),form_low,form_high,gbm.atom.nall());
  } else
    kernel_pack_nbor<numtyp><<<GX,BX,0,gbm.pair_stream>>>
          ((vec4 *)gbm.atom.dev_x.begin(),gbm.nbor.dev_nbor.begin(), 
           gbm.atom.inum(), start, inum,
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
EXTERN void gb_gpu_name(const int id, const int max_nbors, char * name) {
  string sname=GBMF[0].gpu.name(id)+", "+
              gb_gpu_toa(GBMF[0].gpu.cores(id))+" cores, "+
              gb_gpu_toa(GBMF[0].gpu.gigabytes(id))+" GB, "+
              gb_gpu_toa(GBMF[0].gpu.clock_rate(id))+" GHZ";
  strcpy(name,sname.c_str());
}

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
EXTERN bool gb_gpu_init(int &ij_size, const int ntypes, const double gamma,
                  const double upsilon, const double mu, double **shape,
                  double **well, double **cutsq, double **sigma, 
                  double **epsilon, double *host_lshape, int **form,
                  double **host_lj1, double **host_lj2, double **host_lj3, 
                  double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, 
                  const int max_nbors, const int thread, const int gpu_id) {
  assert(thread<MAX_GPU_THREADS);
  
  GBMF[thread].gpu.init();

  if (GBMF[thread].gpu.num_devices()==0)
    return false;                   

  ij_size=IJ_SIZE;
  return GBMF[thread].init(ij_size, ntypes, gamma, upsilon, mu, shape,
                           well, cutsq, sigma, epsilon, host_lshape, form,
                           host_lj1, host_lj2, host_lj3, host_lj4, offset,
                           special_lj, nlocal, nall, max_nbors, false, 
                           gpu_id);
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
EXTERN void gb_gpu_clear(const int thread) {
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
  atom.add_x_data(host_x,host_type);
  atom.add_q_data(host_quat[0]);

  atom.copy_x_data(stream);
  atom.copy_q_data(stream);
  atom.time_atom.stop();
}

EXTERN void gb_gpu_atom(double **host_x, double **host_quat, 
                        const int *host_type, const bool rebuild, 
                        const int thread) {
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
  success=true;
    
  gbm.nbor.time_nbor.start();

  int mn=0;
  for (int i=0; i<inum; i++)
    mn=std::max(mn,numj[i]);

  if (nall>gbm.max_atoms)
    gbm.resize_atom(nall,success);  
  if (nlocal>gbm.max_local || mn>gbm._max_nbors)
    gbm.resize_local(nlocal,mn,success);
  if (!success)
    return false;
    
  gbm.atom.nall(nall);
  gbm.atom.inum(inum);

  if (gbm.multiple_forms) {
    int ij_size=gbm.nbor.host_ij.numel();
    if (inum*2<ij_size) {
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
      gbm.nbor.host_ij.copy_to_device(gbm.nbor.dev_nbor.begin()+inum,
                                      2*inum,gbm.pair_stream);
    } else {
      int p=0, acc=0;
      int offset=0;
      int half=ij_size/2;
      int hi=0;
      for (int i=0; i<inum; i++) {
        int itype=type[ilist[i]];
        if (gbm.host_form[itype][itype]==ELLIPSE_ELLIPSE) {
          gbm.host_olist[p]=ilist[i];
          gbm.nbor.host_ij[hi]=numj[ilist[i]];
          gbm.nbor.host_ij[hi+half]=acc;
          acc+=numj[ilist[i]];
          p++;
          hi++;
          if (hi==half) {
            gbm.nbor.host_ij.copy_to_device(gbm.nbor.dev_nbor.begin()+inum+offset,
                                            half,gbm.pair_stream);
            gbm.nbor.host_ij.copy_to_device(half,gbm.nbor.dev_nbor.begin()+
                                                 inum*2+offset,
                                            half,gbm.pair_stream);
            hi=0;
            offset+=half;
            CUDA_SAFE_CALL(cudaStreamSynchronize(gbm.pair_stream));
          }
        }
      }
      gbm.last_ellipse=p;
      for (int i=0; i<inum; i++) {
        int itype=type[ilist[i]];
        if (gbm.host_form[itype][itype]!=ELLIPSE_ELLIPSE) {
          gbm.host_olist[p]=ilist[i];
          gbm.nbor.host_ij[hi]=numj[ilist[i]];
          gbm.nbor.host_ij[hi+half]=acc;
          acc+=numj[ilist[i]];
          p++;
          hi++;
          if (hi==half) {
            gbm.nbor.host_ij.copy_to_device(gbm.nbor.dev_nbor.begin()+inum+offset,
                                            half,gbm.pair_stream);
            gbm.nbor.host_ij.copy_to_device(half,gbm.nbor.dev_nbor.begin()+
                                                 inum*2+offset,
                                            half,gbm.pair_stream);
            hi=0;
            offset+=half;
            CUDA_SAFE_CALL(cudaStreamSynchronize(gbm.pair_stream));
          }
        }
      }
      gbm.nbor.dev_nbor.copy_from_host(gbm.host_olist.begin(),inum);
      if (hi>0) {
        gbm.nbor.host_ij.copy_to_device(gbm.nbor.dev_nbor.begin()+inum+offset,
                                        hi,gbm.pair_stream);
        gbm.nbor.host_ij.copy_to_device(half,gbm.nbor.dev_nbor.begin()+
                                             inum*2+offset,
                                        hi,gbm.pair_stream);
      }
      gbm.nbor.ij_total=0;
    }
  } else {
    gbm.nbor.reset(inum,ilist,numj,gbm.pair_stream);
    gbm.last_ellipse=inum;
  }

  gbm.nbor.time_nbor.stop();
  
  if (gbm.multiple_forms)
    return gbm.host_olist.begin();
  return ilist;
}

EXTERN int * gb_gpu_reset_nbors(const int nall, const int nlocal, const int inum, 
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
void _gb_gpu_nbors(gbmtyp &gbm, const int *ij, const int num_ij, 
        const bool eflag) {
  gbm.nbor.time_nbor.add_to_total();
  // CUDA_SAFE_CALL(cudaStreamSynchronize(gbm.pair_stream)); // Not if timed
  
  memcpy(gbm.nbor.host_ij.begin(),ij,num_ij*sizeof(int));
  gbm.nbor.time_nbor.start();
  gbm.nbor.add(num_ij,gbm.pair_stream);
  gbm.nbor.time_nbor.stop();
}

EXTERN void gb_gpu_nbors(const int *ij, const int num_ij, const bool eflag,
            const int thread) {
  _gb_gpu_nbors(GBMF[thread],ij,num_ij,eflag);
}


template<class numtyp, class acctyp>
void _gb_gpu_enqueue(GBMT &gbm, const bool eflag, const bool vflag) {
  gbm.atom.time_answer.start();
  gbm.atom.copy_answers(eflag,vflag,gbm.pair_stream);
  gbm.atom.time_answer.stop();
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques for all ij interactions
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void _gb_gpu_gayberne(GBMT &gbm, const bool eflag, const bool vflag, 
                      const bool rebuild) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=BLOCK_1D;
  int ans_pitch=6;
  if (eflag)
    ans_pitch++;
  if (vflag)
    ans_pitch+=6;
  
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
           ((vec4*)gbm.atom.dev_x.begin(), (vec4*)gbm.atom.dev_q.begin(), 
            gbm.gamma_upsilon_mu.begin(), gbm.special_lj.begin(), 
            gbm.nbor.dev_nbor.begin(), gbm.atom.inum(),
            gbm.atom.ans.begin(), ans_pitch,gbm.dev_error.begin(),
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
           ((vec4*)gbm.atom.dev_x.begin(), (vec4*)gbm.atom.dev_q.begin(), 
            gbm.gamma_upsilon_mu.begin(), gbm.special_lj.begin(), 
            gbm.nbor.dev_nbor.begin(), gbm.atom.inum(),
            gbm.atom.ans.begin(), ans_pitch,gbm.dev_error.begin(),
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
           ((vec4*)gbm.atom.dev_x.begin(), gbm.special_lj.begin(), 
            gbm.nbor.dev_nbor.begin(), gbm.atom.inum(), gbm.nbor.ij.begin(),
            gbm.atom.ans.begin(), ans_pitch, gbm.dev_error.begin(), eflag, 
            vflag, gbm.last_ellipse, gbm.atom.inum(), gbm.atom.nall());
      else
        kernel_lj<numtyp,acctyp><<<GX,BX,0,gbm.pair_stream>>>
           ((vec4*)gbm.atom.dev_x.begin(), gbm.special_lj.begin(), 
            gbm.nbor.dev_nbor.begin(), gbm.atom.inum(), gbm.nbor.ij.begin(),
            gbm.atom.ans.begin(), ans_pitch,gbm.dev_error.begin(),
            eflag, vflag, gbm.last_ellipse, gbm.atom.inum(), gbm.atom.nall());
    }
    gbm.time_pair.stop();
  } else {
    gbm.time_kernel.start();
    pack_nbors(gbm, GX, BX, 0, gbm.atom.inum(),SPHERE_SPHERE,ELLIPSE_ELLIPSE);
    gbm.time_kernel.stop();
  
    gbm.time_gayberne.start(); 
    kernel_gayberne<numtyp,acctyp><<<GX,BX,0,gbm.pair_stream>>>
         ((vec4*)gbm.atom.dev_x.begin(), (vec4*)gbm.atom.dev_q.begin(), 
          gbm.gamma_upsilon_mu.begin(), gbm.special_lj.begin(), 
          gbm.nbor.dev_nbor.begin(), gbm.atom.inum(),
          gbm.atom.ans.begin(), ans_pitch, gbm.dev_error.begin(), 
          eflag, vflag, gbm.atom.inum(), gbm.atom.nall());
    gbm.time_gayberne.stop();
  }
}

EXTERN void gb_gpu_gayberne(const bool eflag, const bool vflag, const bool rebuild, 
                            const int thread) {
  _gb_gpu_gayberne<PRECISION,ACC_PRECISION>(GBMF[thread],eflag,vflag,rebuild);
  _gb_gpu_enqueue<PRECISION,ACC_PRECISION>(GBMF[thread],eflag,vflag);
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

  gbm.atom.time_atom.add_to_total();
  gbm.nbor.time_nbor.add_to_total();
  gbm.time_kernel.add_to_total();
  gbm.time_gayberne.add_to_total();
  if (gbm.multiple_forms) {
    gbm.time_kernel2.add_to_total();
    gbm.time_gayberne2.add_to_total();
    gbm.time_pair.add_to_total();
  }      
  CUDA_SAFE_CALL(cudaStreamSynchronize(gbm.pair_stream));
  if (gbm.last_ellipse>gbm.atom.inum()) {
    if (eflag || vflag)
      evdw=gbm.atom.energy_virial(ilist,eflag_atom,vflag_atom,eatom,vatom,virial,
                                  f,tor,gbm.atom.inum());
    else
      gbm.atom.copy_asphere(ilist,f,tor,gbm.atom.inum());
  } else {
    if (eflag || vflag)
      evdw=gbm.atom.energy_virial(ilist,eflag_atom,vflag_atom,eatom,vatom,virial,
                                  f,tor,gbm.last_ellipse);
    else
      gbm.atom.copy_asphere(ilist,f,tor,gbm.last_ellipse);
  }
  gbm.atom.time_answer.add_to_total();
  return evdw;
}

EXTERN double gb_gpu_forces(double **f, double **tor, const int *ilist,
                     const bool eflag, const bool vflag, const bool eflag_atom,
                     const bool vflag_atom, double *eatom, double **vatom,
                     double *virial, const int thread) {
  return _gb_gpu_forces<PRECISION,ACC_PRECISION>
                       (GBMF[thread],f,tor,ilist,eflag,vflag,eflag_atom,
                        vflag_atom,eatom,vatom,virial);
}

EXTERN void gb_gpu_time(const int i) {
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

EXTERN int gb_gpu_num_devices() {
  return GBMF[0].gpu.num_devices();
}

EXTERN double gb_gpu_bytes() {
  return GBMF[0].host_memory_usage();
}

