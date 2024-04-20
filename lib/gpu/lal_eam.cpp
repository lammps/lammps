/***************************************************************************
                                   eam.cpp
                             -------------------
                   Trung Dac Nguyen, W. Michael Brown (ORNL)

  Class for acceleration of the eam pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov nguyentd@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "eam_cl.h"
#elif defined(USE_CUDART)
const char *eam=0;
#else
#include "eam_cubin.h"
#endif

#include "lal_eam.h"
#include <cassert>
namespace LAMMPS_AL {
#define EAMT EAM<numtyp, acctyp>


#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
EAMT::EAM() : BaseAtomic<numtyp,acctyp>(),
  _compiled_energy(false), _allocated(false) {
}

template <class numtyp, class acctyp>
EAMT::~EAM() {
  clear();
}

template <class numtyp, class acctyp>
int EAMT::init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
               int **host_type2z2r, int *host_type2frho,
               double ***host_rhor_spline, double ***host_z2r_spline,
               double ***host_frho_spline, double** host_cutsq, double rdr, double rdrho,
               double rhomax, int nrhor, int nrho, int nz2r, int nfrho, int nr,
               const int nlocal, const int nall, const int max_nbors,
               const int maxspecial, const double cell_size,
               const double gpu_split, FILE *_screen)
{
  int max_shared_types=this->device->max_shared_types();

  int onetype=0;
  #ifdef USE_OPENCL
  for (int i=1; i<ntypes; i++)
    if (host_type2frho[i]>=0 && host_type2frho[i]<=nfrho-1) {
      if (onetype>0)
        onetype=-1;
      else if (onetype==0)
        onetype=i*max_shared_types+i;
    }
  if (onetype<0) onetype=0;
  #endif

  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            gpu_split,_screen,eam,"k_eam",onetype);

  if (success!=0)
    return success;

  // allocate fp

  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;

  _max_fp_size=static_cast<int>(static_cast<double>(ef_nall)*1.10);
  _fp.alloc(_max_fp_size,*(this->ucl_device),UCL_READ_WRITE,UCL_READ_WRITE);

  k_energy.set_function(*(this->pair_program),"k_energy");
  k_energy_fast.set_function(*(this->pair_program),"k_energy_fast");
  fp_tex.get_texture(*(this->pair_program),"fp_tex");
  fp_tex.bind_float(_fp,1);

  #if defined(LAL_OCL_EV_JIT)
  k_energy_fast_noev.set_function(*(this->pair_program_noev),"k_energy_fast");
  #else
  k_energy_sel = &k_energy_fast;
  #endif

  _compiled_energy = true;

  // Initialize timers for selected GPU
  time_pair2.init(*(this->ucl_device));
  time_pair2.zero();

  time_fp1.init(*(this->ucl_device));
  time_fp1.zero();

  time_fp2.init(*(this->ucl_device));
  time_fp2.zero();

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;

  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }

  _ntypes=lj_types;
  _cutforcesq=host_cutforcesq;
  _rdr=rdr;
  _rdrho = rdrho;
  _rhomax=rhomax;
  _nrhor=nrhor;
  _nrho=nrho;
  _nz2r=nz2r;
  _nfrho=nfrho;
  _nr=nr;

  UCL_H_Vec<int2> dview_type(lj_types*lj_types,*(this->ucl_device),
                             UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++) {
    dview_type[i].x=0; dview_type[i].y=0;
  }

  // pack type2rhor and type2z2r
  type2rhor_z2r.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<ntypes; i++) {
    for (int j=0; j<ntypes; j++) {
      dview_type[i*lj_types+j].x=host_type2rhor[i][j];
      dview_type[i*lj_types+j].y=host_type2z2r[i][j];
    }
  }

  ucl_copy(type2rhor_z2r,dview_type,false);

  // pack type2frho
  UCL_H_Vec<int> dview_type2frho(lj_types,*(this->ucl_device),
                                 UCL_WRITE_ONLY);

  type2frho.alloc(lj_types,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<ntypes; i++)
    dview_type2frho[i]=host_type2frho[i];
  ucl_copy(type2frho,dview_type2frho,false);

  // pack frho_spline
  UCL_H_Vec<numtyp4> dview_frho_spline(nfrho*(nrho+1),*(this->ucl_device),
                                       UCL_WRITE_ONLY);

  for (int ix=0; ix<nfrho; ix++)
    for (int iy=0; iy<nrho+1; iy++) {
      dview_frho_spline[ix*(nrho+1)+iy].x=host_frho_spline[ix][iy][0];
      dview_frho_spline[ix*(nrho+1)+iy].y=host_frho_spline[ix][iy][1];
      dview_frho_spline[ix*(nrho+1)+iy].z=host_frho_spline[ix][iy][2];
      dview_frho_spline[ix*(nrho+1)+iy].w=0;
    }

  frho_spline1.alloc(nfrho*(nrho+1),*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(frho_spline1,dview_frho_spline,false);
  frho_spline1_tex.get_texture(*(this->pair_program),"frho_sp1_tex");
  frho_spline1_tex.bind_float(frho_spline1,4);

  for (int ix=0; ix<nfrho; ix++)
    for (int iy=0; iy<nrho+1; iy++) {
      dview_frho_spline[ix*(nrho+1)+iy].x=host_frho_spline[ix][iy][3];
      dview_frho_spline[ix*(nrho+1)+iy].y=host_frho_spline[ix][iy][4];
      dview_frho_spline[ix*(nrho+1)+iy].z=host_frho_spline[ix][iy][5];
      dview_frho_spline[ix*(nrho+1)+iy].w=host_frho_spline[ix][iy][6];
    }

  frho_spline2.alloc(nfrho*(nrho+1),*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(frho_spline2,dview_frho_spline,false);
  frho_spline2_tex.get_texture(*(this->pair_program),"frho_sp2_tex");
  frho_spline2_tex.bind_float(frho_spline2,4);

  // pack rhor_spline
  UCL_H_Vec<numtyp4> dview_rhor_spline(nrhor*(nr+1),*(this->ucl_device),
                                       UCL_WRITE_ONLY);

  for (int ix=0; ix<nrhor; ix++)
    for (int iy=0; iy<nr+1; iy++) {
      dview_rhor_spline[ix*(nr+1)+iy].x=host_rhor_spline[ix][iy][0];
      dview_rhor_spline[ix*(nr+1)+iy].y=host_rhor_spline[ix][iy][1];
      dview_rhor_spline[ix*(nr+1)+iy].z=host_rhor_spline[ix][iy][2];
      dview_rhor_spline[ix*(nr+1)+iy].w=(numtyp)0;
    }

  rhor_spline1.alloc(nrhor*(nr+1),*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(rhor_spline1,dview_rhor_spline,false);
  rhor_spline1_tex.get_texture(*(this->pair_program),"rhor_sp1_tex");
  rhor_spline1_tex.bind_float(rhor_spline1,4);

  for (int ix=0; ix<nrhor; ix++)
    for (int iy=0; iy<nr+1; iy++) {
      dview_rhor_spline[ix*(nr+1)+iy].x=host_rhor_spline[ix][iy][3];
      dview_rhor_spline[ix*(nr+1)+iy].y=host_rhor_spline[ix][iy][4];
      dview_rhor_spline[ix*(nr+1)+iy].z=host_rhor_spline[ix][iy][5];
      dview_rhor_spline[ix*(nr+1)+iy].w=host_rhor_spline[ix][iy][6];
    }

  rhor_spline2.alloc(nrhor*(nr+1),*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(rhor_spline2,dview_rhor_spline,false);
  rhor_spline2_tex.get_texture(*(this->pair_program),"rhor_sp2_tex");
  rhor_spline2_tex.bind_float(rhor_spline2,4);

  // pack z2r_spline
  UCL_H_Vec<numtyp4> dview_z2r_spline(nz2r*(nr+1),*(this->ucl_device),
                                      UCL_WRITE_ONLY);

  for (int ix=0; ix<nz2r; ix++)
    for (int iy=0; iy<nr+1; iy++) {
      dview_z2r_spline[ix*(nr+1)+iy].x=host_z2r_spline[ix][iy][0];
      dview_z2r_spline[ix*(nr+1)+iy].y=host_z2r_spline[ix][iy][1];
      dview_z2r_spline[ix*(nr+1)+iy].z=host_z2r_spline[ix][iy][2];
      dview_z2r_spline[ix*(nr+1)+iy].w=(numtyp)0;
    }

  z2r_spline1.alloc(nz2r*(nr+1),*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(z2r_spline1,dview_z2r_spline,false);
  z2r_spline1_tex.get_texture(*(this->pair_program),"z2r_sp1_tex");
  z2r_spline1_tex.bind_float(z2r_spline1,4);

  for (int ix=0; ix<nz2r; ix++)
    for (int iy=0; iy<nr+1; iy++) {
      dview_z2r_spline[ix*(nr+1)+iy].x=host_z2r_spline[ix][iy][3];
      dview_z2r_spline[ix*(nr+1)+iy].y=host_z2r_spline[ix][iy][4];
      dview_z2r_spline[ix*(nr+1)+iy].z=host_z2r_spline[ix][iy][5];
      dview_z2r_spline[ix*(nr+1)+iy].w=host_z2r_spline[ix][iy][6];
    }

  z2r_spline2.alloc(nz2r*(nr+1),*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(z2r_spline2,dview_z2r_spline,false);
  z2r_spline2_tex.get_texture(*(this->pair_program),"z2r_sp2_tex");
  z2r_spline2_tex.bind_float(z2r_spline2,4);

  UCL_H_Vec<numtyp> host_write(lj_types*lj_types,*(this->ucl_device),
                               UCL_WRITE_ONLY);
  host_write.zero();
  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,cutsq,host_write,host_cutsq);

  _allocated=true;
  this->_max_bytes=type2rhor_z2r.row_bytes()
        + type2frho.row_bytes()
        + rhor_spline1.row_bytes()
        + rhor_spline2.row_bytes()
        + frho_spline1.row_bytes()
        + frho_spline2.row_bytes()
        + z2r_spline1.row_bytes()
        + z2r_spline2.row_bytes()
        + cutsq.row_bytes()
        + _fp.device.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void EAMT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  type2rhor_z2r.clear();
  type2frho.clear();
  rhor_spline1.clear();
  rhor_spline2.clear();
  frho_spline1.clear();
  frho_spline2.clear();
  z2r_spline1.clear();
  z2r_spline2.clear();
  cutsq.clear();

  _fp.clear();

  time_pair2.clear();
  time_fp1.clear();
  time_fp2.clear();

  if (_compiled_energy) {
    k_energy_fast.clear();
    k_energy.clear();
    #if defined(LAL_OCL_EV_JIT)
    k_energy_fast_noev.clear();
    #endif
    _compiled_energy=false;
  }

  this->clear_atomic();
}

template <class numtyp, class acctyp>
double EAMT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(EAM<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then compute per-atom fp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::compute(const int f_ago, const int inum_full, const int nlocal,
                   const int nall, double **host_x, int *host_type,
                   int *ilist, int *numj, int **firstneigh,
                   const bool eflag_in, const bool vflag_in,
                   const bool /*eatom*/, const bool /*vatom*/,
                   int &host_start, const double cpu_time,
                   bool &success, void **fp_ptr) {
  this->acc_timers();
  int eflag, vflag;
  if (eflag_in) eflag=2;
  else eflag=0;
  if (vflag_in) vflag=2;
  else vflag=0;

  this->set_kernel(eflag,vflag);

  if (this->device->time_device()) {
    // Put time from the second part to the total time_pair
    this->time_pair.add_time_to_total(time_pair2.time());

    // Add transfer time from device -> host after part 1
    this->atom->add_transfer_time(time_fp1.time());

    // Add transfer time from host -> device before part 2
    this->atom->add_transfer_time(time_fp2.time());
  }

  // ------------------- Resize FP Array for EAM --------------------

  if (nall>_max_fp_size) {
    _max_fp_size=static_cast<int>(static_cast<double>(nall)*1.10);
    _fp.resize(_max_fp_size);
    fp_tex.bind_float(_fp,1);
  }
  *fp_ptr=_fp.host.begin();

  // ----------------------------------------------------------------

  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    this->resize_atom(0,nall,success);
    this->zero_timers();
    return;
  }

  int ago=this->hd_balancer.ago_first(f_ago);
  int inum=this->hd_balancer.balance(ago,inum_full,cpu_time);
  this->ans->inum(inum);
  host_start=inum;

  // -----------------------------------------------------------------

  if (ago==0) {
    this->reset_nbors(nall, inum, ilist, numj, firstneigh, success);
    if (!success)
      return;
  }

  this->atom->cast_x_data(host_x,host_type);
  this->atom->add_x_data(host_x,host_type);

  loop(eflag,vflag);

  // copy fp from device to host for comm
  _nlocal=nlocal;
  time_fp1.start();
  _fp.update_host(nlocal,true);
  time_fp1.stop();
  time_fp1.sync_stop();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU and then compute per-atom fp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** EAMT::compute(const int ago, const int inum_full, const int nall,
                    double **host_x, int *host_type, double *sublo,
                    double *subhi, tagint *tag, int **nspecial,
                    tagint **special, const bool eflag_in,
                    const bool vflag_in, const bool /*eatom*/,
                    const bool /*vatom*/, int &host_start, int **ilist, int **jnum,
                    const double cpu_time, bool &success, int &inum,
                    void **fp_ptr) {
  this->acc_timers();
  int eflag, vflag;
  if (eflag_in) eflag=2;
  else eflag=0;
  if (vflag_in) vflag=2;
  else vflag=0;

  this->set_kernel(eflag,vflag);

  if (this->device->time_device()) {
    // Put time from the second part to the total time_pair
    this->time_pair.add_time_to_total(time_pair2.time());

    // Add transfer time from device -> host after part 1
    this->atom->add_transfer_time(time_fp1.time());

    // Add transfer time from host -> device before part 2
    this->atom->add_transfer_time(time_fp2.time());
  }

  // ------------------- Resize FP Array for EAM --------------------

  if (nall>_max_fp_size) {
    _max_fp_size=static_cast<int>(static_cast<double>(nall)*1.10);
    _fp.resize(_max_fp_size);
    fp_tex.bind_float(_fp,1);
  }
  *fp_ptr=_fp.host.begin();

  // -----------------------------------------------------------------

  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    this->resize_atom(0,nall,success);
    this->zero_timers();
    return nullptr;
  }

  // load balance, returning the atom count on the device (inum)
  this->hd_balancer.balance(cpu_time);
  inum=this->hd_balancer.get_gpu_count(ago,inum_full);
  this->ans->inum(inum);
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    this->build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                          sublo, subhi, tag, nspecial, special, success);
    if (!success)
      return nullptr;
  } else {
    this->atom->cast_x_data(host_x,host_type);
    this->atom->add_x_data(host_x,host_type);
  }
  *ilist=this->nbor->host_ilist.begin();
  *jnum=this->nbor->host_acc.begin();

  loop(eflag,vflag);

  // copy fp from device to host for comm
  _nlocal=inum_full;
  time_fp1.start();
  _fp.update_host(inum_full,true);
  time_fp1.stop();
  time_fp1.sync_stop();

  return this->nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Update per-atom fp, and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::compute2(int *ilist, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom) {
  if (this->ans->inum()==0)
    return;

  this->hd_balancer.start_timer();
  time_fp2.start();
  this->add_fp_data();
  time_fp2.stop();

  loop2(eflag,vflag);
  if (ilist == nullptr)
    this->ans->copy_answers(eflag,vflag,eatom,vatom, this->ans->inum());
  else
    this->ans->copy_answers(eflag,vflag,eatom,vatom, ilist, this->ans->inum());

  this->device->add_ans_object(this->ans);
  this->hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Calculate per-atom fp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int EAMT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();

  if (shared_types) {
    #if defined(LAL_OCL_EV_JIT)
    if (eflag || vflag) k_energy_sel = &k_energy_fast;
    else k_energy_sel = &k_energy_fast_noev;
    #endif

    k_energy_sel->set_size(GX,BX);
    k_energy_sel->run(&this->atom->x, &type2rhor_z2r, &type2frho,
                      &rhor_spline2, &frho_spline1, &frho_spline2, &cutsq,
                      &this->nbor->dev_nbor,  &this->_nbor_data->begin(),
                      &_fp, &this->ans->engv, &eflag, &ainum,
                      &nbor_pitch, &_ntypes, &_cutforcesq, &_rdr, &_rdrho,
                      &_rhomax, &_nrho, &_nr, &this->_threads_per_atom);
  } else {
    this->k_energy.set_size(GX,BX);
    this->k_energy.run(&this->atom->x, &type2rhor_z2r, &type2frho,
                       &rhor_spline2, &frho_spline1, &frho_spline2, &cutsq,
                       &this->nbor->dev_nbor, &this->_nbor_data->begin(), &_fp,
                       &this->ans->engv,&eflag, &ainum, &nbor_pitch,
                       &_ntypes, &_cutforcesq, &_rdr, &_rdrho, &_rhomax, &_nrho,
                       &_nr, &this->_threads_per_atom);
  }

  this->time_pair.stop();
  return ainum;
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::loop2(const bool _eflag, const bool _vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;

  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair2.start();

  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &_fp, &type2rhor_z2r,
                          &rhor_spline1, &z2r_spline1, &z2r_spline2, &cutsq,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &_cutforcesq, &_rdr,
                          &_nr, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &_fp, &type2rhor_z2r, &rhor_spline1,
                     &z2r_spline1, &z2r_spline2, &cutsq, &this->nbor->dev_nbor,
                     &this->_nbor_data->begin(), &this->ans->force,
                     &this->ans->engv, &eflag, &vflag, &ainum, &nbor_pitch,
                     &_ntypes, &_cutforcesq, &_rdr, &_nr,
                     &this->_threads_per_atom);
  }

  this->time_pair2.stop();
}

template class EAM<PRECISION,ACC_PRECISION>;
}
