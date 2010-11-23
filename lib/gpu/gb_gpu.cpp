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

#include <iostream>
#include <cassert>
#include <math.h>

#include "gb_gpu_memory.h"

using namespace std;

static GB_GPU_Memory<PRECISION,ACC_PRECISION> GBMF;
#define GBMT GB_GPU_Memory<numtyp,acctyp>

template<class numtyp, class acctyp>
void gb_gpu_pack_nbors(GBMT &gbm, const int GX, const int BX, const int start, 
                const int inum, const int form_low, const int form_high) {
  int stride=gbm.nbor->nbor_pitch();
  int anall=gbm.atom->nall();
  if (gbm.shared_types) {
    GBMF.k_gb_nbor_fast.set_size(GX,BX);
    GBMF.k_gb_nbor_fast.run(&gbm.atom->dev_x.begin(),
              &gbm.cut_form.begin(), &gbm.nbor->dev_nbor.begin(), &stride,
              &start, &inum, &gbm.nbor->dev_packed.begin(), &form_low,
              &form_high, &anall);
  } else {
    GBMF.k_gb_nbor.set_size(GX,BX);
    GBMF.k_gb_nbor.run(&gbm.atom->dev_x.begin(), &gbm.cut_form.begin(),
              &gbm._lj_types, &gbm.nbor->dev_nbor.begin(), &stride,
              &start, &inum, &gbm.nbor->dev_packed.begin(), &form_low,
              &form_high, &anall);
  }
}

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
bool gb_gpu_init(const int ntypes, const double gamma,
                 const double upsilon, const double mu, double **shape,
                 double **well, double **cutsq, double **sigma,
                 double **epsilon, double *host_lshape, int **form,
                 double **host_lj1, double **host_lj2, double **host_lj3,
                 double **host_lj4, double **offset, double *special_lj,
                 const int inum, const int nall, const int max_nbors, 
                 const double cell_size, int &gpu_mode, FILE *screen) {
  GBMF.clear();
  gpu_mode=GBMF.device->gpu_mode();
  double gpu_split=GBMF.device->particle_split();
  int first_gpu=GBMF.device->first_device();
  int last_gpu=GBMF.device->last_device();
  int world_me=GBMF.device->world_me();
  int gpu_rank=GBMF.device->gpu_rank();
  int procs_per_gpu=GBMF.device->procs_per_gpu();

  GBMF.device->init_message(screen,"gayberne",first_gpu,last_gpu);

  bool message=false;
  if (world_me==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  if (world_me==0) {
    bool init_ok=GBMF.init(ntypes, gamma, upsilon, mu, shape, well, cutsq, 
                           sigma, epsilon, host_lshape, form, host_lj1, 
                           host_lj2, host_lj3, host_lj4, offset, special_lj, 
                           inum, nall, max_nbors, cell_size, gpu_split, screen);
    if (!init_ok)
      return false;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (message)
    fprintf(screen,"Done.\n");
        
  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing GPU %d on core %d...",gpu_rank,i);
      else
        fprintf(screen,"Initializing GPUs %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0) {
      bool init_ok=GBMF.init(ntypes, gamma, upsilon, mu, shape, well, cutsq, 
                             sigma, epsilon, host_lshape, form, host_lj1, 
                             host_lj2, host_lj3, host_lj4, offset, special_lj, 
                             inum, nall, max_nbors, cell_size, gpu_split, 
                             screen);
      if (!init_ok)
        return false;
    }
    MPI_Barrier(GBMF.device->gpu_comm);
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");
  return true;
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
void gb_gpu_clear() {
  GBMF.clear();
}

// ---------------------------------------------------------------------------
// Build neighbor list on device
// ---------------------------------------------------------------------------
template <class gbmtyp>
inline void _gb_gpu_build_nbor_list(gbmtyp &gbm, const int inum,
                                    const int host_inum, const int nall, 
                                    double **host_x, double **host_quat,
                                    int *host_type, double *boxlo,
                                    double *boxhi, bool &success) {
  gbm.nbor_time_avail=true;

  success=true;
  gbm.resize_atom(inum,nall,success);
  gbm.resize_local(inum,host_inum,gbm.nbor->max_nbors(),0,success);
  if (!success)
    return;
    
  gbm.atom->cast_copy_x(host_x,host_type);
  int mn;
  gbm.nbor->build_nbor_list(inum, host_inum, nall, *gbm.atom,
                            boxlo, boxhi, NULL, NULL, NULL, success, mn);
  gbm.nbor->copy_unpacked(inum,mn);
  gbm.last_ellipse=inum;
  gbm.max_last_ellipse=inum;
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host and (if spheres) reorder so ellipses first
// ---------------------------------------------------------------------------
template <class gbmtyp>
void _gb_gpu_reset_nbors(gbmtyp &gbm, const int nall,
                          const int inum, const int osize,
                          int *ilist, int *numj,
                          int *type, int **firstneigh,
                          bool &success) {
  success=true;
    
  gbm.nbor_time_avail=true;

  int mn=gbm.nbor->max_nbor_loop(inum,numj);
  gbm.resize_atom(inum,nall,success);
  gbm.resize_local(inum,0,mn,osize,success);
  if (!success)
    return;
    
  if (gbm.multiple_forms) {
    int p=0;
    for (int i=0; i<osize; i++) {
      int itype=type[ilist[i]];
      if (gbm.host_form[itype][itype]==ELLIPSE_ELLIPSE) {
        gbm.host_olist[p]=ilist[i];
        p++;
      }
    }
    gbm.max_last_ellipse=p;
    gbm.last_ellipse=std::min(inum,gbm.max_last_ellipse);
    for (int i=0; i<osize; i++) {
      int itype=type[ilist[i]];
      if (gbm.host_form[itype][itype]!=ELLIPSE_ELLIPSE) {
        gbm.host_olist[p]=ilist[i];
        p++;
      }
    }
    gbm.nbor->get_host(inum,gbm.host_olist.begin(),numj,firstneigh,
                      gbm.block_size());
    gbm.nbor->copy_unpacked(inum,mn);
    return;
  }
  gbm.last_ellipse=inum;
  gbm.max_last_ellipse=inum;
  gbm.nbor->get_host(inum,ilist,numj,firstneigh,gbm.block_size());
  gbm.nbor->copy_unpacked(inum,mn);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void _gb_gpu_gayberne(GBMT &gbm, const bool _eflag, const bool _vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=gbm.block_size();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;
  
  int GX=static_cast<int>(ceil(static_cast<double>(gbm.atom->inum())/BX));
  int stride=gbm.nbor->nbor_pitch();
  int ainum=gbm.atom->inum();
  int anall=gbm.atom->nall();

  if (gbm.multiple_forms) {
    gbm.time_kernel.start();
    if (gbm.last_ellipse>0) {
      // ------------ ELLIPSE_ELLIPSE and ELLIPSE_SPHERE ---------------
      GX=static_cast<int>(ceil(static_cast<double>(gbm.last_ellipse)/
                               static_cast<double>(BX)));
      gb_gpu_pack_nbors(gbm,GX,BX, 0, gbm.last_ellipse,ELLIPSE_SPHERE,
			ELLIPSE_ELLIPSE);
      gbm.time_kernel.stop();

      gbm.time_gayberne.start();
      GBMF.k_gayberne.set_size(GX,BX);
      GBMF.k_gayberne.run(&gbm.atom->dev_x.begin(),
           &gbm.atom->dev_quat.begin(), &gbm.shape.begin(), &gbm.well.begin(),
           &gbm.gamma_upsilon_mu.begin(), &gbm.sigma_epsilon.begin(), 
           &gbm._lj_types, &gbm.lshape.begin(), &gbm.nbor->dev_nbor.begin(),
           &stride, &gbm.atom->dev_ans.begin(),&ainum,&gbm.atom->dev_engv.begin(),
           &gbm.dev_error.begin(), &eflag, &vflag, &gbm.last_ellipse, &anall);
      gbm.time_gayberne.stop();

      if (gbm.last_ellipse==gbm.atom->inum()) {
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
      GX=static_cast<int>(ceil(static_cast<double>(gbm.atom->inum()-
                               gbm.last_ellipse)/BX));
      gb_gpu_pack_nbors(gbm,GX,BX,gbm.last_ellipse,gbm.atom->inum(),
			SPHERE_ELLIPSE,SPHERE_ELLIPSE);
      gbm.time_kernel2.stop();

      gbm.time_gayberne2.start();
      GBMF.k_sphere_gb.set_size(GX,BX);
      GBMF.k_sphere_gb.run(&gbm.atom->dev_x.begin(),&gbm.atom->dev_quat.begin(),
              &gbm.shape.begin(), &gbm.well.begin(), 
              &gbm.gamma_upsilon_mu.begin(), &gbm.sigma_epsilon.begin(), 
              &gbm._lj_types, &gbm.lshape.begin(), 
              &gbm.nbor->dev_nbor.begin(), &stride, &gbm.atom->dev_ans.begin(),
              &gbm.atom->dev_engv.begin(), &gbm.dev_error.begin(), &eflag,
              &vflag, &gbm.last_ellipse, &ainum, &anall);
      gbm.time_gayberne2.stop();
   } else {
      gbm.atom->dev_ans.zero();
      gbm.atom->dev_engv.zero();
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
    if (gbm.last_ellipse<gbm.atom->inum()) {
      if (gbm.shared_types) {
        GBMF.k_lj_fast.set_size(GX,BX);
        GBMF.k_lj_fast.run(&gbm.atom->dev_x.begin(), &gbm.lj1.begin(),
                           &gbm.lj3.begin(), &gbm.gamma_upsilon_mu.begin(),
                           &stride, &gbm.nbor->dev_packed.begin(),
                           &gbm.atom->dev_ans.begin(),
                           &gbm.atom->dev_engv.begin(), &gbm.dev_error.begin(),
                           &eflag, &vflag, &gbm.last_ellipse, &ainum, &anall);
      } else {
        GBMF.k_lj.set_size(GX,BX);
        GBMF.k_lj.run(&gbm.atom->dev_x.begin(), &gbm.lj1.begin(),
                      &gbm.lj3.begin(), &gbm._lj_types, 
                      &gbm.gamma_upsilon_mu.begin(), &stride, 
                      &gbm.nbor->dev_packed.begin(), &gbm.atom->dev_ans.begin(),
                      &gbm.atom->dev_engv.begin(), &gbm.dev_error.begin(),
                      &eflag, &vflag, &gbm.last_ellipse, &ainum, &anall);
      }
    }
    gbm.time_pair.stop();
  } else {
    gbm.time_kernel.start();
    gb_gpu_pack_nbors(gbm, GX, BX, 0, gbm.atom->inum(),SPHERE_SPHERE,
		      ELLIPSE_ELLIPSE);
    gbm.time_kernel.stop();
    gbm.time_gayberne.start(); 
    GBMF.k_gayberne.set_size(GX,BX);
    GBMF.k_gayberne.run(&gbm.atom->dev_x.begin(), &gbm.atom->dev_quat.begin(),
            &gbm.shape.begin(), &gbm.well.begin(), 
            &gbm.gamma_upsilon_mu.begin(), &gbm.sigma_epsilon.begin(), 
            &gbm._lj_types, &gbm.lshape.begin(), &gbm.nbor->dev_nbor.begin(),
            &stride, &gbm.atom->dev_ans.begin(), &ainum,
            &gbm.atom->dev_engv.begin(), &gbm.dev_error.begin(),
            &eflag, &vflag, &ainum, &anall);
    gbm.time_gayberne.stop();
  }
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, torques, energies
// ---------------------------------------------------------------------------
template <class gbmtyp>
inline int * _gb_gpu_compute_n(gbmtyp &gbm, const int timestep, const int ago,
		               const int inum_full, const int nall,
			       double **host_x, int *host_type,
			       double *boxlo, double *boxhi, const bool eflag,
			       const bool vflag, const bool eatom,
                               const bool vatom, int &host_start,
		               const double cpu_time, bool &success,
			       double **host_quat) {
  gbm.acc_timers();
  if (inum_full==0) {
    gbm.zero_timers();
    return NULL;
  }

  gbm.hd_balancer.balance(cpu_time,gbm.nbor->gpu_nbor());
  int inum=gbm.hd_balancer.get_gpu_count(timestep,ago,inum_full);
  gbm.atom->inum(inum);
  gbm.last_ellipse=std::min(inum,gbm.max_last_ellipse);
  host_start=inum;
  
  // Build neighbor list on GPU if necessary
  if (ago==0) {
    _gb_gpu_build_nbor_list(gbm, inum, inum_full-inum, nall, host_x,
                            host_quat, host_type, boxlo, boxhi, success);
    if (!success)
      return NULL;
    gbm.atom->cast_quat_data(host_quat[0]);
    gbm.hd_balancer.start_timer();
  } else {    
    gbm.atom->cast_x_data(host_x,host_type);
    gbm.atom->cast_quat_data(host_quat[0]);
    gbm.hd_balancer.start_timer();
    gbm.atom->add_x_data(host_x,host_type);
  }

  gbm.atom->add_other_data();

  _gb_gpu_gayberne<PRECISION,ACC_PRECISION>(gbm,eflag,vflag);
  gbm.atom->copy_answers(eflag,vflag,eatom,vatom);
  gbm.hd_balancer.stop_timer();
  return gbm.device->nbor.host_nbor.begin();
}

int * gb_gpu_compute_n(const int timestep, const int ago, const int inum_full,
	 	       const int nall, double **host_x, int *host_type,
                       double *boxlo, double *boxhi, const bool eflag,
		       const bool vflag, const bool eatom, const bool vatom,
                       int &host_start, const double cpu_time, bool &success,
		       double **host_quat) {
  return _gb_gpu_compute_n(GBMF, timestep, ago, inum_full, nall, host_x,
			   host_type, boxlo, boxhi, eflag, vflag, eatom, vatom,
                           host_start, cpu_time, success, host_quat);
}  

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, torques,..
// ---------------------------------------------------------------------------
template <class gbmtyp>
inline int * _gb_gpu_compute(gbmtyp &gbm, const int timestep, const int f_ago,
			     const int inum_full,const int nall,double **host_x,
			     int *host_type, int *ilist, int *numj,
			     int **firstneigh, const bool eflag,
			     const bool vflag, const bool eatom,
                             const bool vatom, int &host_start,
			     const double cpu_time, bool &success,
			     double **host_quat) {
  gbm.acc_timers();
  if (inum_full==0) {
    gbm.zero_timers();
    return NULL;
  }
  
  int ago=gbm.hd_balancer.ago_first(f_ago);
  int inum=gbm.hd_balancer.balance(timestep,ago,inum_full,cpu_time,
				   gbm.nbor->gpu_nbor());
  gbm.atom->inum(inum);
  gbm.last_ellipse=std::min(inum,gbm.max_last_ellipse);
  host_start=inum;

  if (ago==0) {
    _gb_gpu_reset_nbors(gbm, nall, inum, inum_full, ilist, numj, host_type,
		        firstneigh, success);
    if (!success)
      return NULL;
  }
  int *list;
  if (gbm.multiple_forms)
    list=gbm.host_olist.begin();
  else
    list=ilist;

  gbm.atom->cast_x_data(host_x,host_type);
  gbm.atom->cast_quat_data(host_quat[0]);
  gbm.hd_balancer.start_timer();
  gbm.atom->add_x_data(host_x,host_type);
  gbm.atom->add_other_data();

  _gb_gpu_gayberne<PRECISION,ACC_PRECISION>(gbm,eflag,vflag);
  gbm.atom->copy_answers(eflag,vflag,eatom,vatom,list);
  gbm.hd_balancer.stop_timer();
  return list;
}

int * gb_gpu_compute(const int timestep, const int ago, const int inum_full,
	 	     const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh,
		     const bool eflag, const bool vflag, const bool eatom,
                     const bool vatom, int &host_start, const double cpu_time,
                     bool &success, double **host_quat) {
  return _gb_gpu_compute(GBMF, timestep, ago, inum_full, nall, host_x,
			 host_type, ilist, numj, firstneigh, eflag, vflag,
			 eatom, vatom, host_start, cpu_time, success,
                         host_quat);
}

// ---------------------------------------------------------------------------
// Return memory usage
// ---------------------------------------------------------------------------
double gb_gpu_bytes() {
  return GBMF.host_memory_usage();
}
