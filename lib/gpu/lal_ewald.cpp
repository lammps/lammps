/***************************************************************************
                                lal_ewald.cpp
                             -------------------
                            W. Michael Brown (ORNL)
                            Axel Kohlmeyer (Temple)

  Class for Ewald (reciprocal-space) acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : developers@lammps.org
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "ewald_cl.h"
#elif defined(USE_CUDART)
const char *ewald=0;
#else
#include "ewald_cubin.h"
#endif

#include "lal_ewald.h"
#include <cassert>
#include <cmath>

namespace LAMMPS_AL {
#define EwaldGPUT EwaldGPU<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
EwaldGPUT::EwaldGPU() : _allocated(false), _compiled(false), _max_bytes(0) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  ewald_program=nullptr;
  ucl_device=nullptr;
  _kmax=0;
  _kcount=0;
  _cs_kmax=-1;
  _cs_nlocal=0;
}

template <class numtyp, class acctyp>
EwaldGPUT::~EwaldGPU() {
  clear(0.0);
  k_cssn.clear();
  k_structure.clear();
  k_field.clear();
  if (ewald_program) delete ewald_program;
  delete ans;
}

template <class numtyp, class acctyp>
int EwaldGPUT::init(const int nlocal, const int nall, FILE *_screen,
                    int &flag) {
  _max_bytes=10;
  screen=_screen;

  flag=device->init(*ans,nlocal,nall);
  if (flag!=0)
    return flag;
  if (device->ptx_arch()>0.0 && device->ptx_arch()<1.1) {
    flag=-4;
    return flag;
  }

  if (ucl_device!=device->gpu) _compiled=false;
  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pair_block_size();

  compile_kernels(*ucl_device);

  time_in.init(*ucl_device);
  time_in.zero();
  time_out.init(*ucl_device);
  time_out.zero();
  time_map.init(*ucl_device);
  time_map.zero();
  time_rho.init(*ucl_device);
  time_rho.zero();
  time_interp.init(*ucl_device);
  time_interp.zero();
  _cpu_idle_time=0.0;

  pos_tex.bind_float(atom->x,4);
  q_tex.bind_float(atom->q,1);

  _allocated=true;
  _max_bytes=0;
  _max_an_bytes=ans->gpu_bytes();
  _cs_kmax=-1;
  _cs_nlocal=0;
  _kmax=0;
  _kcount=0;
  return flag;
}

// ---------------------------------------------------------------------------
// Upload the (constant per box) k-vectors and grid parameters
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EwaldGPUT::setup(const int kmax, const int kcount, int *kxvecs,
                      int *kyvecs, int *kzvecs, double *ug, double **eg,
                      double **vg, double *unitk, bool &success) {
  _kmax=kmax;
  _kcount=kcount;
  _unitk[0]=(numtyp)unitk[0];
  _unitk[1]=(numtyp)unitk[1];
  _unitk[2]=(numtyp)unitk[2];

  // (re)allocate and upload the integer k-vectors (length kcount)

  d_kxvecs.clear();
  d_kyvecs.clear();
  d_kzvecs.clear();
  success = success && (d_kxvecs.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);
  success = success && (d_kyvecs.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);
  success = success && (d_kzvecs.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);

  // (re)allocate and upload the energy/field/virial coefficients

  d_eg.clear();
  d_ug.clear();
  d_vg.clear();
  success = success && (d_eg.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);
  success = success && (d_ug.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);
  success = success && (d_vg.alloc(kcount*6,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);

  // (re)allocate the structure-factor buffers (host+device, length kcount)

  d_sfacrl.clear();
  d_sfacim.clear();
  success = success && (d_sfacrl.alloc(kcount,*ucl_device,UCL_READ_WRITE,
                                       UCL_READ_WRITE)==UCL_SUCCESS);
  success = success && (d_sfacim.alloc(kcount,*ucl_device,UCL_READ_WRITE,
                                       UCL_READ_WRITE)==UCL_SUCCESS);
  if (!success)
    return;

  UCL_H_Vec<int> view;
  view.view(kxvecs,kcount,*ucl_device);
  ucl_copy(d_kxvecs,view,false);
  view.view(kyvecs,kcount,*ucl_device);
  ucl_copy(d_kyvecs,view,false);
  view.view(kzvecs,kcount,*ucl_device);
  ucl_copy(d_kzvecs,view,false);

  UCL_H_Vec<numtyp4> eg_view;
  eg_view.alloc(kcount,*ucl_device);
  UCL_H_Vec<numtyp> ug_view, vg_view;
  ug_view.alloc(kcount,*ucl_device);
  vg_view.alloc(kcount*6,*ucl_device);
  for (int k=0; k<kcount; k++) {
    eg_view[k].x=(numtyp)eg[k][0];
    eg_view[k].y=(numtyp)eg[k][1];
    eg_view[k].z=(numtyp)eg[k][2];
    eg_view[k].w=(numtyp)0.0;
    ug_view[k]=(numtyp)ug[k];
    for (int j=0; j<6; j++) vg_view[k*6+j]=(numtyp)vg[k][j];
  }
  ucl_copy(d_eg,eg_view,false);
  ucl_copy(d_ug,ug_view,false);
  ucl_copy(d_vg,vg_view,false);
  eg_view.clear();
  ug_view.clear();
  vg_view.clear();

  // force reallocation of the per-atom cs/sn arrays (kmax may have changed)

  _cs_kmax=-1;
  _cs_nlocal=0;
}

// ---------------------------------------------------------------------------
// K-space contribution to the per-atom field and force
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EwaldGPUT::compute_forces(double *host_sfacrl_all, double *host_sfacim_all,
                               const double qscale, const int slabflag,
                               const int eflag_atom, const int vflag_atom,
                               double *host_eatom, double **host_vatom,
                               bool &success) {
  if (_nlocal==0)
    return;

  // upload the global structure factors (after the host MPI_Allreduce)
  for (int k=0; k<_kcount; k++) {
    d_sfacrl[k]=(acctyp)host_sfacrl_all[k];
    d_sfacim[k]=(acctyp)host_sfacim_all[k];
  }
  d_sfacrl.update_device(_kcount,false);
  d_sfacim.update_device(_kcount,false);

  _qscale=(numtyp)qscale;
  _slabflag=slabflag;
  _eflag_atom=eflag_atom;
  _vflag_atom=vflag_atom;

  time_interp.start();
  const int BX=_block_size;
  const int GX=static_cast<int>(ceil(static_cast<double>(_nlocal)/BX));
  k_field.set_size(GX,BX);
  k_field.run(&atom->q, &d_cs, &d_sn, &d_kxvecs, &d_kyvecs, &d_kzvecs, &d_eg,
              &d_ug, &d_vg, &d_sfacrl, &d_sfacim, &ans->force, &d_eatom,
              &d_vatom, &_qscale, &_slabflag, &_eflag_atom, &_vflag_atom,
              &_kmax, &_nlocal, &_kcount);
  time_interp.stop();

  ans->copy_answers(false,false,false,false,0);
  device->add_ans_object(ans);

  // copy the raw per-atom energy/virial back to the host
  if (eflag_atom) {
    d_eatom.update_host(_nlocal,false);
    for (int i=0; i<_nlocal; i++)
      host_eatom[i]=d_eatom[i];
  }
  if (vflag_atom) {
    d_vatom.update_host(_nlocal*6,false);
    for (int i=0; i<_nlocal; i++)
      for (int j=0; j<6; j++)
        host_vatom[i][j]=d_vatom[i*6+j];
  }

  success=true;
}

// ---------------------------------------------------------------------------
// (Re)allocate the per-atom cos/sin arrays for the current kmax and nlocal
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EwaldGPUT::resize_cssn(const int nlocal, bool &success) {
  if (nlocal<=_cs_nlocal && _kmax==_cs_kmax)
    return;

  int newn=static_cast<int>(static_cast<double>(nlocal)*1.10)+1;
  int sz=(2*_kmax+1)*3*newn;

  d_cs.clear();
  d_sn.clear();
  d_eatom.clear();
  d_vatom.clear();
  success = success && (d_cs.alloc(sz,*ucl_device)==UCL_SUCCESS);
  success = success && (d_sn.alloc(sz,*ucl_device)==UCL_SUCCESS);
  success = success && (d_eatom.alloc(newn,*ucl_device,UCL_WRITE_ONLY,
                                      UCL_READ_WRITE)==UCL_SUCCESS);
  success = success && (d_vatom.alloc(newn*6,*ucl_device,UCL_WRITE_ONLY,
                                      UCL_READ_WRITE)==UCL_SUCCESS);
  if (!success)
    return;

  _cs_nlocal=newn;
  _cs_kmax=_kmax;
}

// ---------------------------------------------------------------------------
// Compute the local (per-rank) structure factors on the device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int EwaldGPUT::structure(const int ago, const int nlocal, const int nall,
                         double **host_x, int *host_type, double *host_q,
                         double *host_sfacrl, double *host_sfacim,
                         bool &success) {
  // roll up the previous step's per-phase timers into their totals
  acc_timers();

  _nlocal=nlocal;

  // a rank with no local atoms contributes zero to the structure factors;
  // zero the host arrays so the subsequent MPI_Allreduce is correct
  if (nlocal==0) {
    for (int k=0; k<_kcount; k++) {
      host_sfacrl[k]=0.0;
      host_sfacim[k]=0.0;
    }
    return 0;
  }

  ans->inum(nlocal);

  if (ago==0) {
    resize_atom(nlocal,nall,success);
    if (!success)
      return 0;
  }
  resize_cssn(nlocal,success);
  if (!success)
    return 0;

  time_in.start();
  atom->cast_x_data(host_x,host_type);
  atom->cast_q_data(host_q);
  atom->add_x_data(host_x,host_type);
  atom->add_q_data();
  time_in.stop();

  const int BX=_block_size;

  // cos/sin recurrence: one thread per atom
  time_map.start();
  const int GX=static_cast<int>(ceil(static_cast<double>(nlocal)/BX));
  k_cssn.set_size(GX,BX);
  k_cssn.run(&atom->x, &d_cs, &d_sn, &_unitk[0], &_unitk[1], &_unitk[2],
             &_kmax, &_nlocal);
  time_map.stop();

  // structure factor reduction: one block per k-vector
  time_rho.start();
  k_structure.set_size(_kcount,BX);
  k_structure.run(&atom->q, &d_cs, &d_sn, &d_kxvecs, &d_kyvecs, &d_kzvecs,
                  &d_sfacrl, &d_sfacim, &_kmax, &_nlocal, &_kcount);
  time_rho.stop();

  time_out.start();
  d_sfacrl.update_host(_kcount,false);
  d_sfacim.update_host(_kcount,false);
  time_out.stop();
  for (int k=0; k<_kcount; k++) {
    host_sfacrl[k]=d_sfacrl[k];
    host_sfacim[k]=d_sfacim[k];
  }

  return 0;
}

template <class numtyp, class acctyp>
void EwaldGPUT::clear(const double cpu_time) {
  if (!_allocated)
    return;
  _allocated=false;

  d_kxvecs.clear();
  d_kyvecs.clear();
  d_kzvecs.clear();
  d_eg.clear();
  d_ug.clear();
  d_vg.clear();
  d_cs.clear();
  d_sn.clear();
  d_eatom.clear();
  d_vatom.clear();
  d_sfacrl.clear();
  d_sfacim.clear();

  acc_timers();
  device->output_kspace_times(time_in,time_out,time_map,time_rho,time_interp,
                              *ans,_max_bytes+_max_an_bytes,cpu_time,
                              _cpu_idle_time,screen);

  time_in.clear();
  time_out.clear();
  time_map.clear();
  time_rho.clear();
  time_interp.clear();

  ans->clear();
}

template <class numtyp, class acctyp>
double EwaldGPUT::host_memory_usage() const {
  return device->atom.host_memory_usage()+
         sizeof(EwaldGPU<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void EwaldGPUT::compile_kernels(UCL_Device &dev) {
  if (_compiled)
    return;

  std::string flags=device->compile_string();

  if (ewald_program) delete ewald_program;
  ewald_program=new UCL_Program(dev);
  ewald_program->load_string(ewald,flags.c_str(),nullptr,screen);

  k_cssn.set_function(*ewald_program,"k_ewald_cssn");
  k_structure.set_function(*ewald_program,"k_ewald_structure");
  k_field.set_function(*ewald_program,"k_ewald_field");
  pos_tex.get_texture(*ewald_program,"pos_tex");
  q_tex.get_texture(*ewald_program,"q_tex");

  _compiled=true;
}

template class EwaldGPU<PRECISION,ACC_PRECISION>;
}
