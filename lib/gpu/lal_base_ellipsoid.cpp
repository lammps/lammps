/***************************************************************************
                              base_ellipsoid.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Base class for acceleration of ellipsoid potentials

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Thu May 5 2011
    email                : brownw@ornl.gov
 ***************************************************************************/

#include "lal_base_ellipsoid.h"
#include <cstdlib>
using namespace LAMMPS_AL;

#if defined(USE_OPENCL)
#include "ellipsoid_nbor_cl.h"
#elif defined(USE_CUDART)
const char *ellipsoid_nbor=0;
#else
#include "ellipsoid_nbor_cubin.h"
#endif

#define BaseEllipsoidT BaseEllipsoid<numtyp, acctyp>
extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
BaseEllipsoidT::BaseEllipsoid() : _compiled(false), _max_bytes(0) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  nbor=new Neighbor();
}

template <class numtyp, class acctyp>
BaseEllipsoidT::~BaseEllipsoid() {
  delete ans;
  delete nbor;
}

template <class numtyp, class acctyp>
int BaseEllipsoidT::bytes_per_atom(const int max_nbors) const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+
         nbor->bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
int BaseEllipsoidT::init_base(const int nlocal, const int nall,
                              const int max_nbors, const int maxspecial,
                              const double cell_size, const double gpu_split,
                              FILE *_screen, const int ntypes, int **h_form,
                              const void *ellipsoid_program,
                              const void *lj_program, const char *k_name,
                              const bool ellip_sphere) {
  screen=_screen;
  _ellipsoid_sphere=ellip_sphere;

  int gpu_nbor=0;
  if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=1;
  else if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_HYB_NEIGH)
    gpu_nbor=2;

  int _gpu_host=0;
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_split,gpu_nbor);
  if (host_nlocal>0)
    _gpu_host=1;

  _threads_per_atom=device->threads_per_atom();

  int success=device->init(*ans,false,true,nlocal,host_nlocal,nall,nbor,
                           maxspecial,_gpu_host,max_nbors,cell_size,true,
                           1);
  if (success!=0)
    return success;

  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->block_ellipse();
  compile_kernels(*ucl_device,ellipsoid_program,lj_program,k_name,ellip_sphere);

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_nbor,gpu_split);

  // Initialize timers for the selected GPU
  time_lj.init(*ucl_device);
  time_nbor1.init(*ucl_device);
  time_ellipsoid.init(*ucl_device);
  time_nbor2.init(*ucl_device);
  time_ellipsoid2.init(*ucl_device);
  time_nbor3.init(*ucl_device);
  time_ellipsoid3.init(*ucl_device);
  time_lj.zero();
  time_nbor1.zero();
  time_ellipsoid.zero();
  time_nbor2.zero();
  time_ellipsoid2.zero();
  time_nbor3.zero();
  time_ellipsoid3.zero();

  // See if we want fast GB-sphere or sphere-sphere calculations
  _host_form=h_form;
  _multiple_forms=false;
  for (int i=1; i<ntypes; i++)
    for (int j=i; j<ntypes; j++)
      if (_host_form[i][j]!=ELLIPSE_ELLIPSE)
        _multiple_forms=true;
  if (_multiple_forms && host_nlocal>0)
    return -8;
  if (_multiple_forms && gpu_nbor!=0)
    return -9;

  if (_multiple_forms)
    ans->force.zero();

  // Memory for ilist ordered by particle type
  if (host_olist.alloc(nbor->max_atoms(),*ucl_device)!=UCL_SUCCESS)
    return -3;

  _max_an_bytes=ans->gpu_bytes()+nbor->gpu_bytes();

  neigh_tex.bind_float(atom->x,4);
  pos_tex.bind_float(atom->x,4);
  quat_tex.bind_float(atom->quat,4);
  lj_pos_tex.bind_float(atom->x,4);
  lj_quat_tex.bind_float(atom->quat,4);

  return 0;
}

template <class numtyp, class acctyp>
void BaseEllipsoidT::estimate_gpu_overhead() {
  device->estimate_gpu_overhead(2,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void BaseEllipsoidT::clear_base() {
  // Output any timing information
  output_times();
  host_olist.clear();

  if (_compiled) {
    k_nbor_fast.clear();
    k_nbor.clear();
    k_ellipsoid.clear();
    k_ellipsoid_sphere.clear();
    k_sphere_ellipsoid.clear();
    k_lj_fast.clear();
    k_lj.clear();
    delete nbor_program;
    delete ellipsoid_program;
    delete lj_program;
    _compiled=false;
  }

  time_nbor1.clear();
  time_ellipsoid.clear();
  time_nbor2.clear();
  time_ellipsoid2.clear();
  time_nbor3.clear();
  time_ellipsoid3.clear();
  time_lj.clear();
  hd_balancer.clear();

  nbor->clear();
  ans->clear();
  device->clear();
}

template <class numtyp, class acctyp>
void BaseEllipsoidT::output_times() {
  // Output any timing information
  acc_timers();
  double single[10], times[10];

  single[0]=atom->transfer_time()+ans->transfer_time();
  single[1]=nbor->time_nbor.total_seconds()+nbor->time_hybrid1.total_seconds()+
            nbor->time_hybrid2.total_seconds();
  single[2]=time_nbor1.total_seconds()+time_nbor2.total_seconds()+
            time_nbor3.total_seconds()+nbor->time_nbor.total_seconds();
  single[3]=time_ellipsoid.total_seconds()+time_ellipsoid2.total_seconds()+
            time_ellipsoid3.total_seconds();
  if (_multiple_forms)
    single[4]=time_lj.total_seconds();
  else
    single[4]=0;
  single[5]=atom->cast_time()+ans->cast_time();
  single[6]=_gpu_overhead;
  single[7]=_driver_overhead;
  single[8]=ans->cpu_idle_time();
  single[9]=nbor->bin_time();

  MPI_Reduce(single,times,10,MPI_DOUBLE,MPI_SUM,0,device->replica());
  double avg_split=hd_balancer.all_avg_split();

  _max_bytes+=atom->max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,
             device->replica());
  double max_mb=mpi_max_bytes/(1024*1024);
  double t_time=times[0]+times[1]+times[2]+times[3]+times[4]+times[5];

  if (device->replica_me()==0)
    if (screen && times[5]>0.0) {
      int replica_size=device->replica_size();

      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"    Device Time Info (average): ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (device->procs_per_gpu()==1 && t_time>0) {
        fprintf(screen,"Data Transfer:   %.4f s.\n",times[0]/replica_size);
        fprintf(screen,"Data Cast/Pack:  %.4f s.\n",times[5]/replica_size);
        fprintf(screen,"Neighbor copy:   %.4f s.\n",times[1]/replica_size);
        if (nbor->gpu_nbor()>0)
          fprintf(screen,"Neighbor build:  %.4f s.\n",times[2]/replica_size);
        else
          fprintf(screen,"Neighbor unpack: %.4f s.\n",times[2]/replica_size);
        fprintf(screen,"Force calc:      %.4f s.\n",times[3]/replica_size);
        fprintf(screen,"LJ calc:         %.4f s.\n",times[4]/replica_size);
      }
      if (nbor->gpu_nbor()==2)
        fprintf(screen,"Neighbor (CPU):  %.4f s.\n",times[9]/replica_size);
      if (times[6]>0)
        fprintf(screen,"Device Overhead: %.4f s.\n",times[6]/replica_size);
      fprintf(screen,"Average split:   %.4f.\n",avg_split);
      fprintf(screen,"Threads / atom:  %d.\n",_threads_per_atom);
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);
      fprintf(screen,"CPU Driver_Time: %.4f s.\n",times[7]/replica_size);
      fprintf(screen,"CPU Idle_Time:   %.4f s.\n",times[8]/replica_size);
      fprintf(screen,"-------------------------------------");
      fprintf(screen,"--------------------------------\n\n");
    }
  _max_bytes=0.0;
}

// ---------------------------------------------------------------------------
// Pack neighbors to limit thread divergence for lj-lj and ellipse
// ---------------------------------------------------------------------------
template<class numtyp, class acctyp>
void BaseEllipsoidT::pack_nbors(const int GX, const int BX, const int start,
                                const int inum, const int form_low,
                                const int form_high, const bool shared_types,
                                int ntypes) {
  int stride=nbor->nbor_pitch();
  if (shared_types) {
    k_nbor_fast.set_size(GX,BX);
    k_nbor_fast.run(&atom->x, &cut_form, &nbor->dev_nbor, &stride, &start,
                    &inum, &nbor->dev_packed, &form_low, &form_high);
  } else {
    k_nbor.set_size(GX,BX);
    k_nbor.run(&atom->x, &cut_form, &ntypes, &nbor->dev_nbor, &stride,
               &start, &inum, &nbor->dev_packed, &form_low, &form_high);
  }
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseEllipsoidT::reset_nbors(const int nall, const int inum,
                                 const int osize, int *ilist,
                                 int *numj, int *type, int **firstneigh,
                                 bool &success) {
  success=true;

  int mn=nbor->max_nbor_loop(osize,numj,ilist);
  resize_atom(nall,success);
  resize_local(inum,0,mn,osize,success);
  if (!success)
    return;

  if (_multiple_forms) {
    int p=0;
    for (int i=0; i<osize; i++) {
      int itype=type[ilist[i]];
      if (_host_form[itype][itype]==ELLIPSE_ELLIPSE) {
        host_olist[p]=ilist[i];
        p++;
      }
    }
    _max_last_ellipse=p;
    _last_ellipse=std::min(inum,_max_last_ellipse);
    for (int i=0; i<osize; i++) {
      int itype=type[ilist[i]];
      if (_host_form[itype][itype]!=ELLIPSE_ELLIPSE) {
        host_olist[p]=ilist[i];
        p++;
      }
    }
    nbor->get_host(inum,host_olist.begin(),numj,firstneigh,block_size());
    nbor->copy_unpacked(inum,mn);
    return;
  }
  _last_ellipse=inum;
  _max_last_ellipse=inum;
  nbor->get_host(inum,ilist,numj,firstneigh,block_size());
  nbor->copy_unpacked(inum,mn);

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
}

// ---------------------------------------------------------------------------
// Build neighbor list on device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline void BaseEllipsoidT::build_nbor_list(const int inum, const int host_inum,
                                            const int nall, double **host_x,
                                            int *host_type, double *sublo,
                                            double *subhi, tagint *tag,
                                            int **nspecial, tagint **special,
                                            bool &success) {
  success=true;
  resize_atom(nall,success);
  resize_local(inum,host_inum,nbor->max_nbors(),0,success);
  if (!success)
    return;
  atom->cast_copy_x(host_x,host_type);

  int mn;
  nbor->build_nbor_list(host_x, inum, host_inum, nall, *atom, sublo, subhi, tag,
                        nspecial, special, success, mn);
  nbor->copy_unpacked(inum,mn);
  _last_ellipse=inum;
  _max_last_ellipse=inum;

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int* BaseEllipsoidT::compute(const int f_ago, const int inum_full,
                             const int nall, double **host_x, int *host_type,
                             int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag,
                             const bool eatom, const bool vatom,
                             int &host_start, const double cpu_time,
                             bool &success, double **host_quat) {
  acc_timers();
  if (inum_full==0) {
    host_start=0;
    zero_timers();
    return NULL;
  }

  int ago=hd_balancer.ago_first(f_ago);
  int inum=hd_balancer.balance(ago,inum_full,cpu_time);
  ans->inum(inum);
  _last_ellipse=std::min(inum,_max_last_ellipse);
  host_start=inum;

  if (ago==0) {
    reset_nbors(nall, inum, inum_full, ilist, numj, host_type, firstneigh,
                success);
    if (!success)
      return NULL;
  }
  int *list;
  if (_multiple_forms)
    list=host_olist.begin();
  else
    list=ilist;

  atom->cast_x_data(host_x,host_type);
  atom->cast_quat_data(host_quat[0]);
  hd_balancer.start_timer();
  atom->add_x_data(host_x,host_type);
  atom->add_quat_data();

  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom,list);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();
  return list;
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** BaseEllipsoidT::compute(const int ago, const int inum_full, const int nall,
                              double **host_x, int *host_type, double *sublo,
                              double *subhi, tagint *tag, int **nspecial,
                              tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom,
                              int &host_start, int **ilist, int **jnum,
                              const double cpu_time, bool &success,
                              double **host_quat) {
  acc_timers();
  if (inum_full==0) {
    host_start=0;
    zero_timers();
    return NULL;
  }

  hd_balancer.balance(cpu_time);
  int inum=hd_balancer.get_gpu_count(ago,inum_full);
  ans->inum(inum);
  _last_ellipse=std::min(inum,_max_last_ellipse);
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, success);
    if (!success)
      return NULL;
    atom->cast_quat_data(host_quat[0]);
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    atom->cast_quat_data(host_quat[0]);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }

  atom->add_quat_data();
  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();
  return nbor->host_jlist.begin()-host_start;
}

template <class numtyp, class acctyp>
double BaseEllipsoidT::host_memory_usage_base() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(BaseEllipsoid<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void BaseEllipsoidT::compile_kernels(UCL_Device &dev,
                                     const void *ellipsoid_string,
                                     const void *lj_string,
                                     const char *kname, const bool e_s) {
  if (_compiled)
    return;

  std::string kns=kname;
  std::string s_sphere_ellipsoid=kns+"_sphere_ellipsoid";
  std::string s_ellipsoid_sphere=kns+"_ellipsoid_sphere";
  std::string s_lj=kns+"_lj";
  std::string s_lj_fast=kns+"_lj_fast";

  std::string flags=device->compile_string();

  nbor_program=new UCL_Program(dev);
  nbor_program->load_string(ellipsoid_nbor,flags.c_str());
  k_nbor_fast.set_function(*nbor_program,"kernel_nbor_fast");
  k_nbor.set_function(*nbor_program,"kernel_nbor");
  neigh_tex.get_texture(*nbor_program,"pos_tex");

  ellipsoid_program=new UCL_Program(dev);
  ellipsoid_program->load_string(ellipsoid_string,flags.c_str());
  k_ellipsoid.set_function(*ellipsoid_program,kname);
  pos_tex.get_texture(*ellipsoid_program,"pos_tex");
  quat_tex.get_texture(*ellipsoid_program,"quat_tex");

  lj_program=new UCL_Program(dev);
  lj_program->load_string(lj_string,flags.c_str());
  k_sphere_ellipsoid.set_function(*lj_program,s_sphere_ellipsoid.c_str());
  k_lj_fast.set_function(*lj_program,s_lj_fast.c_str());
  k_lj.set_function(*lj_program,s_lj.c_str());
  if (e_s)
    k_ellipsoid_sphere.set_function(*lj_program,s_ellipsoid_sphere.c_str());
  lj_pos_tex.get_texture(*lj_program,"pos_tex");
  lj_quat_tex.get_texture(*lj_program,"quat_tex");

  _compiled=true;
}

template class BaseEllipsoid<PRECISION,ACC_PRECISION>;

