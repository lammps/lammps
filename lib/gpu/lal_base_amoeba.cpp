/***************************************************************************
                               base_amoeba.cpp
                             -------------------
                            Trung Dac Nguyen (Northwestern)

  Base class for pair styles needing per-particle data for position,
  charge, and type.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#include "lal_base_amoeba.h"

namespace LAMMPS_AL {
#define BaseAmoebaT BaseAmoeba<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
BaseAmoebaT::BaseAmoeba() : _compiled(false), _max_bytes(0), short_nbor_polar_avail(false) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  nbor=new Neighbor();
  pair_program=nullptr;
  ucl_device=nullptr;
}

template <class numtyp, class acctyp>
BaseAmoebaT::~BaseAmoeba() {
  delete ans;
  delete nbor;
  k_multipole.clear();
  k_udirect2b.clear();
  k_umutual2b.clear();
  k_fphi_uind.clear();
  k_fphi_mpole.clear();
  k_polar.clear();
  k_special15.clear();
  k_short_nbor.clear();

  #if 0 // !defined(USE_OPENCL) && !defined(USE_HIP)
  if (fft_plan_created) cufftDestroy(plan);
  #endif

  if (pair_program) delete pair_program;
}

template <class numtyp, class acctyp>
int BaseAmoebaT::bytes_per_atom_atomic(const int max_nbors) const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+
         nbor->bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
int BaseAmoebaT::init_atomic(const int nlocal, const int nall,
                             const int max_nbors, const int maxspecial,
                             const int maxspecial15,
                             const double cell_size, const double gpu_split,
                             FILE *_screen, const void *pair_program,
                             const char *k_name_multipole,
                             const char *k_name_udirect2b,
                             const char *k_name_umutual2b,
                             const char *k_name_polar,
                             const char *k_name_fphi_uind,
                             const char *k_name_fphi_mpole,
                             const char *k_name_short_nbor,
                             const char* k_name_special15) {
  screen=_screen;

  int gpu_nbor=0;
  if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=1;
  else if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_HYB_NEIGH)
    gpu_nbor=2;

  int _gpu_host=0;
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_split,gpu_nbor);
  if (host_nlocal>0)
    _gpu_host=1;

  _threads_per_atom=device->threads_per_charge();

  bool charge = true;
  bool rot = false;
  bool vel = false;
  _extra_fields = 24; // round up to accomodate quadruples of numtyp values
                      // rpole 13; uind 3; uinp 3; amtype, amgroup; pval
  int success=device->init(*ans,charge,rot,nlocal,nall,maxspecial,vel,_extra_fields/4);
  if (success!=0)
    return success;

  if (ucl_device!=device->gpu) _compiled=false;

  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pair_block_size();
  _block_bio_size=device->block_bio_pair();
  compile_kernels(*ucl_device,pair_program,k_name_multipole,
                   k_name_udirect2b, k_name_umutual2b,k_name_polar,
                   k_name_fphi_uind, k_name_fphi_mpole,
                   k_name_short_nbor, k_name_special15);

  if (_threads_per_atom>1 && gpu_nbor==0) {
    nbor->packing(true);
    _nbor_data=&(nbor->dev_packed);
  } else {
    _nbor_data=&(nbor->dev_nbor);
  }

  bool alloc_packed=false;
  success = device->init_nbor(nbor,nlocal,host_nlocal,nall,maxspecial,
                              _gpu_host,max_nbors,cell_size,alloc_packed,
                              _threads_per_atom);
  if (success!=0)
    return success;

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_nbor,gpu_split);

  // Initialize timers for the selected GPU
  time_pair.init(*ucl_device);
  time_pair.zero();

  pos_tex.bind_float(atom->x,4);
  q_tex.bind_float(atom->q,1);

  _max_an_bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  _maxspecial=maxspecial;
  _maxspecial15=maxspecial15;

  // allocate per-atom array tep

  int ef_nall=nlocal; //nall;
  if (ef_nall==0)
    ef_nall=2000;

  dev_short_nbor.alloc(ef_nall*(2+max_nbors),*(this->ucl_device),UCL_READ_WRITE);

  _max_tep_size=static_cast<int>(static_cast<double>(ef_nall)*1.10);
  _tep.alloc(_max_tep_size*4,*(this->ucl_device),UCL_READ_WRITE,UCL_READ_WRITE);

  _max_fieldp_size = _max_tep_size;
  _fieldp.alloc(_max_fieldp_size*8,*(this->ucl_device),UCL_READ_WRITE,UCL_READ_WRITE);

  _max_thetai_size = 0;

  _nmax = nall;
  dev_nspecial15.alloc(nall,*(this->ucl_device),UCL_READ_ONLY);
  dev_special15.alloc(_maxspecial15*nall,*(this->ucl_device),UCL_READ_ONLY);
  dev_special15_t.alloc(nall*_maxspecial15,*(this->ucl_device),UCL_READ_ONLY);

  #if 0 // !defined(USE_OPENCL) && !defined(USE_HIP)
  fft_plan_created = false;
  #endif

  #ifdef ASYNC_DEVICE_COPY
  _end_command_queue=ucl_device->num_queues();
  ucl_device->push_command_queue();
  #endif

  return success;
}

template <class numtyp, class acctyp>
void BaseAmoebaT::estimate_gpu_overhead(const int add_kernels) {
  device->estimate_gpu_overhead(1+add_kernels,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void BaseAmoebaT::clear_atomic() {
  // Output any timing information
  acc_timers();
  double avg_split=hd_balancer.all_avg_split();
  _gpu_overhead*=hd_balancer.timestep();
  _driver_overhead*=hd_balancer.timestep();
  device->output_times(time_pair,*ans,*nbor,avg_split,_max_bytes+_max_an_bytes,
                       _gpu_overhead,_driver_overhead,_threads_per_atom,screen);

  time_pair.clear();
  hd_balancer.clear();

  dev_short_nbor.clear();
  nbor->clear();
  ans->clear();

  _tep.clear();
  _fieldp.clear();
  _thetai1.clear();
  _thetai2.clear();
  _thetai3.clear();
  _igrid.clear();
  _fdip_phi1.clear();
  _fdip_phi2.clear();
  _fdip_sum_phi.clear();
  _cgrid_brick.clear();

  dev_nspecial15.clear();
  dev_special15.clear();
  dev_special15_t.clear();
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int * BaseAmoebaT::reset_nbors(const int nall, const int inum, int *ilist,
                                   int *numj, int **firstneigh, bool &success) {
  success=true;

  int mn=nbor->max_nbor_loop(inum,numj,ilist);
  resize_atom(inum,nall,success);
  resize_local(inum,mn,success);
  if (!success)
    return nullptr;

  nbor->get_host(inum,ilist,numj,firstneigh,block_size());

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;

  return ilist;
}

// ---------------------------------------------------------------------------
// Build neighbor list on device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline int BaseAmoebaT::build_nbor_list(const int inum, const int host_inum,
                                        const int nall, double **host_x,
                                        int *host_type, double *sublo,
                                        double *subhi, tagint *tag,
                                        int **nspecial, tagint **special,
                                        int *nspecial15, tagint **special15,
                                        bool &success) {
  success=true;
  resize_atom(inum,nall,success);
  resize_local(inum,host_inum,nbor->max_nbors(),success);
  if (!success)
    return 0;
  atom->cast_copy_x(host_x,host_type);

  int mn;
  nbor->build_nbor_list(host_x, inum, host_inum, nall, *atom, sublo, subhi,
                        tag, nspecial, special, success, mn, ans->error_flag);

  // add one-five neighbors

  if (_maxspecial15>0) {
    UCL_H_Vec<int> view_nspecial15;
    UCL_H_Vec<tagint> view_special15;
    view_nspecial15.view(nspecial15,nall,*ucl_device);
    view_special15.view(special15[0],nall*_maxspecial15,*ucl_device);
    ucl_copy(dev_nspecial15,view_nspecial15,nall,false);
    ucl_copy(dev_special15_t,view_special15,_maxspecial15*nall,false);
    nbor->transpose(dev_special15, dev_special15_t, _maxspecial15, nall);

    add_onefive_neighbors();
  }

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
  return mn;
}

// ---------------------------------------------------------------------------
// Prepare for multiple kernel calls in a time step:
//   - reallocate per-atom arrays, if needed
//   - transfer extra data from host to device
//   - build the full neighbor lists for use by different kernels
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
int** BaseAmoebaT::precompute(const int ago, const int inum_full, const int nall,
                              double **host_x, int *host_type, int *host_amtype,
                              int *host_amgroup, double **host_rpole,
                              double **host_uind, double **host_uinp, double *host_pval,
                              double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special,
                              int *nspecial15, tagint **special15,
                              const bool eflag_in, const bool vflag_in,
                              const bool eatom, const bool vatom, int &host_start,
                              int **&ilist, int **&jnum, const double cpu_time,
                              bool &success, double *host_q, double * /*boxlo*/, double * /*prd*/) {
  acc_timers();
  if (eatom) _eflag=2;
  else if (eflag_in) _eflag=1;
  else _eflag=0;
  if (vatom) _vflag=2;
  else if (vflag_in) _vflag=1;
  else _vflag=0;

  #ifdef LAL_NO_BLOCK_REDUCE
  if (_eflag) _eflag=2;
  if (_vflag) _vflag=2;
  #endif

  set_kernel(_eflag,_vflag);

  // ------------------- Resize 1-5 neighbor arrays ------------------------

  if (nall>_nmax) {
    _nmax = nall;
    dev_nspecial15.clear();
    dev_special15.clear();
    dev_special15_t.clear();
    dev_nspecial15.alloc(nall,*(this->ucl_device),UCL_READ_ONLY);
    dev_special15.alloc(_maxspecial15*nall,*(this->ucl_device),UCL_READ_ONLY);
    dev_special15_t.alloc(nall*_maxspecial15,*(this->ucl_device),UCL_READ_ONLY);
  }

  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    resize_atom(0,nall,success);
    zero_timers();
    return nullptr;
  }

  hd_balancer.balance(cpu_time);
  int inum=hd_balancer.get_gpu_count(ago,inum_full);
  ans->inum(inum);
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    _max_nbors = build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, nspecial15, special15,
                    success);
    if (!success)
      return nullptr;
    atom->cast_q_data(host_q);
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    atom->cast_q_data(host_q);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }
  atom->add_q_data();
  cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, host_pval);
  atom->add_extra_data();

  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  // re-allocate dev_short_nbor if necessary
  if (inum_full*(2+_max_nbors) > dev_short_nbor.cols()) {
    int _nmax=static_cast<int>(static_cast<double>(inum_full)*1.10);
    dev_short_nbor.resize((2+_max_nbors)*_nmax);
  }

  hd_balancer.stop_timer();

  return nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Compute multipole real-space part
//   precompute() should be already invoked before mem (re)allocation
//   this is the first part in a time step done on the GPU for AMOEBA for now
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseAmoebaT::compute_multipole_real(const int /*ago*/, const int inum_full,
                                         const int /*nall*/, double ** /*host_x*/,
                                         int * /*host_type*/, int * /*host_amtype*/,
                                         int * /*host_amgroup*/, double ** /*host_rpole*/,
                                         double */*host_pval*/, double * /*sublo*/,
                                         double * /*subhi*/, tagint * /*tag*/,
                                         int ** /*nspecial*/, tagint ** /*special*/,
                                         int * /*nspecial15*/, tagint ** /*special15*/,
                                         const bool /*eflag_in*/, const bool /*vflag_in*/,
                                         const bool /*eatom*/, const bool /*vatom*/,
                                         int & /*host_start*/, int ** /*ilist*/, int ** /*jnum*/,
                                         const double /*cpu_time*/, bool & /*success*/,
                                         const double aewald, const double felec,
                                         const double off2_mpole, double * /*host_q*/,
                                         double * /*boxlo*/, double * /*prd*/, void **tep_ptr) {
  // ------------------- Resize _tep array ------------------------

  if (inum_full>_max_tep_size) {
    _max_tep_size=static_cast<int>(static_cast<double>(inum_full)*1.10);
    _tep.resize(_max_tep_size*4);
  }
  *tep_ptr=_tep.host.begin();

  _off2_mpole = off2_mpole;
  _felec = felec;
  _aewald = aewald;
  multipole_real(_eflag,_vflag);

  // leave the answers (forces, energies and virial) on the device,
  //   only copy them back in the last kernel (polar_real)
  //ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  //device->add_ans_object(ans);

  // copy tep from device to host

  _tep.update_host(_max_tep_size*4,false);
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary, and then compute the direct real space part
//    of the permanent field
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseAmoebaT::compute_udirect2b(int *host_amtype, int *host_amgroup, double **host_rpole,
                                     double **host_uind, double **host_uinp, double *host_pval,
                                     const double aewald, const double off2_polar,
                                     void** fieldp_ptr) {
  // all the necessary data arrays are already copied from host to device

  cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, host_pval);
  atom->add_extra_data();

  *fieldp_ptr=_fieldp.host.begin();

  // specify the correct cutoff and alpha values
  _off2_polar = off2_polar;
  _aewald = aewald;
  udirect2b(_eflag,_vflag);

  // copy field and fieldp from device to host (_fieldp store both arrays, one after another)

  _fieldp.update_host(_max_fieldp_size*8,false);
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary, and then compute the direct real space part
//    of the induced field
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseAmoebaT::compute_umutual2b(int *host_amtype, int *host_amgroup, double ** /*host_rpole*/,
                                    double **host_uind, double **host_uinp, double * /*host_pval*/,
                                    const double aewald, const double off2_polar,
                                    void** /*fieldp_ptr*/) {
  // only copy the necessary data arrays that are updated over the iterations
  // use nullptr for the other arrays that are already copied from host to device
  cast_extra_data(host_amtype, host_amgroup, nullptr, host_uind, host_uinp, nullptr);
  atom->add_extra_data();

  // set the correct cutoff and alpha
  _off2_polar = off2_polar;
  _aewald = aewald;
  // launch the kernel
  umutual2b(_eflag,_vflag);

  // copy field and fieldp from device to host (_fieldp store both arrays, one after another)
  // NOTE: move this step to update_fieldp() to delay device-host transfer
  //       after umutual1 and self are done on the GPU
  // *fieldp_ptr=_fieldp.host.begin();
  // _fieldp.update_host(_max_fieldp_size*8,false);
}

// ---------------------------------------------------------------------------
// Prepare for umutual1() after bspline_fill() is done on host
//   - reallocate per-atom arrays, thetai1, thetai2, thetai3, and igrid if needed
//     host_thetai1, host_thetai2, host_thetai3 are allocated with nmax by bsordermax by 4
//     host_igrid is allocated with nmax by 4
//   - transfer extra data from host to device
// NOTE: can be re-used for fphi_mpole() but with a different bsorder value
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
void BaseAmoebaT::precompute_kspace(const int inum_full, const int bsorder,
                                    double ***host_thetai1, double ***host_thetai2,
                                    double ***host_thetai3, int** host_igrid,
                                    const int nzlo_out, const int nzhi_out,
                                    const int nylo_out, const int nyhi_out,
                                    const int nxlo_out, const int nxhi_out) {
  // update bsorder with that of the kspace solver
  _bsorder = bsorder;

  // allocate or resize per-atom arrays
  // _max_thetai_size, _max_tep_size and _max_fieldp_size are essentially _nmax
  //   will be consolidated once all terms are ready

  if (_max_thetai_size == 0) {
    _max_thetai_size = static_cast<int>(static_cast<double>(inum_full)*1.10);
    _thetai1.alloc(_max_thetai_size*bsorder,*(this->ucl_device),UCL_WRITE_ONLY,UCL_READ_ONLY);
    _thetai2.alloc(_max_thetai_size*bsorder,*(this->ucl_device),UCL_WRITE_ONLY,UCL_READ_ONLY);
    _thetai3.alloc(_max_thetai_size*bsorder,*(this->ucl_device),UCL_WRITE_ONLY,UCL_READ_ONLY);
    _igrid.alloc(_max_thetai_size*4,*(this->ucl_device),UCL_WRITE_ONLY,UCL_READ_ONLY);

    _fdip_phi1.alloc(_max_thetai_size*10,*(this->ucl_device),UCL_READ_WRITE);
    _fdip_phi2.alloc(_max_thetai_size*10,*(this->ucl_device),UCL_READ_WRITE);
    _fdip_sum_phi.alloc(_max_thetai_size*20,*(this->ucl_device),UCL_READ_WRITE);
  } else {
    if ((int)_thetai1.cols()<_max_thetai_size*bsorder) {
      _max_thetai_size=static_cast<int>(static_cast<double>(inum_full)*1.10);
      _thetai1.resize(_max_thetai_size*bsorder);
      _thetai2.resize(_max_thetai_size*bsorder);
      _thetai3.resize(_max_thetai_size*bsorder);
      _igrid.resize(_max_thetai_size*4);

      _fdip_phi1.resize(_max_thetai_size*10);
      _fdip_phi2.resize(_max_thetai_size*10);
      _fdip_sum_phi.resize(_max_thetai_size*20);
    }
  }

  #ifdef ASYNC_DEVICE_COPY
  _thetai1.cq(ucl_device->cq(_end_command_queue));
  _thetai2.cq(ucl_device->cq(_end_command_queue));
  _thetai3.cq(ucl_device->cq(_end_command_queue));
  #endif

  // pack host data to device

  for (int i = 0; i < inum_full; i++)
    for (int j = 0; j < bsorder; j++) {
      int idx = i*bsorder + j;
      numtyp4 v;
      v.x = host_thetai1[i][j][0];
      v.y = host_thetai1[i][j][1];
      v.z = host_thetai1[i][j][2];
      v.w = host_thetai1[i][j][3];
      _thetai1[idx] = v;
    }
  _thetai1.update_device(true);

  for (int i = 0; i < inum_full; i++)
    for (int j = 0; j < bsorder; j++) {
      int idx = i*bsorder + j;
      numtyp4 v;
      v.x = host_thetai2[i][j][0];
      v.y = host_thetai2[i][j][1];
      v.z = host_thetai2[i][j][2];
      v.w = host_thetai2[i][j][3];
      _thetai2[idx] = v;
    }
  _thetai2.update_device(true);

  for (int i = 0; i < inum_full; i++)
    for (int j = 0; j < bsorder; j++) {
      int idx = i*bsorder + j;
      numtyp4 v;
      v.x = host_thetai3[i][j][0];
      v.y = host_thetai3[i][j][1];
      v.z = host_thetai3[i][j][2];
      v.w = host_thetai3[i][j][3];
      _thetai3[idx] = v;
    }
  _thetai3.update_device(true);

  for (int i = 0; i < inum_full; i++) {
    int idx = i*4;
    _igrid[idx+0] = host_igrid[i][0];
    _igrid[idx+1] = host_igrid[i][1];
    _igrid[idx+2] = host_igrid[i][2];
  }
  _igrid.update_device(true);

  // _cgrid_brick holds the grid-based potential

  _nzlo_out = nzlo_out;
  _nzhi_out = nzhi_out;
  _nylo_out = nylo_out;
  _nyhi_out = nyhi_out;
  _nxlo_out = nxlo_out;
  _nxhi_out = nxhi_out;
  _ngridz = nzhi_out - nzlo_out + 1;
  _ngridy = nyhi_out - nylo_out + 1;
  _ngridx = nxhi_out - nxlo_out + 1;
  _num_grid_points = _ngridx * _ngridy * _ngridz;

  int numel = _num_grid_points;
  if (_cgrid_brick.cols() == 0) {
    int nsize=(int)(((double)numel)*1.1);
    _cgrid_brick.alloc(nsize, *(this->ucl_device), UCL_READ_WRITE, UCL_READ_ONLY);
  } else if (numel > (int)_cgrid_brick.cols()) {
    _cgrid_brick.resize(numel);
  }
}

// ---------------------------------------------------------------------------
// fphi_uind = induced potential from grid
// fphi_uind extracts the induced dipole potential from the particle mesh Ewald grid
// NOTE: host_grid_brick is from ic_kspace post_convolution()
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
void BaseAmoebaT::compute_fphi_uind(double ****host_grid_brick,
                                    void **host_fdip_phi1,
                                    void **host_fdip_phi2,
                                    void **host_fdip_sum_phi)
{
  int n = 0;
  for (int iz = _nzlo_out; iz <= _nzhi_out; iz++)
    for (int iy = _nylo_out; iy <= _nyhi_out; iy++)
      for (int ix = _nxlo_out; ix <= _nxhi_out; ix++) {
        numtyp2 v;
        v.x = host_grid_brick[iz][iy][ix][0];
        v.y = host_grid_brick[iz][iy][ix][1];
        _cgrid_brick[n] = v;
        n++;
      }
  _cgrid_brick.update_device(_num_grid_points, false);

  #ifdef ASYNC_DEVICE_COPY
  ucl_device->sync();
  #endif

  // launch the kernel with its execution configuration (see below)
  fphi_uind();

  // copy data from device to host
  _fdip_phi1.update_host(_max_thetai_size*10, false);
  _fdip_phi2.update_host(_max_thetai_size*10, false);
  _fdip_sum_phi.update_host(_max_thetai_size*20, false);

  // return the pointers to the host-side arrays
  *host_fdip_phi1 = _fdip_phi1.host.begin();
  *host_fdip_phi2 = _fdip_phi2.host.begin();
  *host_fdip_sum_phi = _fdip_sum_phi.host.begin();
}

// ---------------------------------------------------------------------------
// Interpolate the potential from the PME grid
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int BaseAmoebaT::fphi_uind() {
  int ainum=ans->inum();
  if (ainum == 0)
    return 0;

  // Compute the block size and grid size to keep all cores busy

  const int BX=block_size();
  const int GX=static_cast<int>(ceil(static_cast<double>(ainum)/BX));

  time_pair.start();
  int ngridxy = _ngridx * _ngridy;
  k_fphi_uind.set_size(GX,BX);
  k_fphi_uind.run(&_thetai1, &_thetai2, &_thetai3, &_igrid, &_cgrid_brick,
                  &_fdip_phi1, &_fdip_phi2, &_fdip_sum_phi, &_bsorder, &ainum,
                  &_nzlo_out, &_nylo_out, &_nxlo_out, &ngridxy, &_ngridx);
  time_pair.stop();

  return GX;
}

// ---------------------------------------------------------------------------
// fphi_mpole = multipole potential from grid (limited to polar_kspace for now)
// fphi_mpole extracts the permanent multipole potential from
//   the particle mesh Ewald grid
// NOTE: host_grid_brick is from p_kspace post_convolution()
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
void BaseAmoebaT::compute_fphi_mpole(double ***host_grid_brick, void **host_fphi, const double felec)
{
  int n = 0;
  for (int iz = _nzlo_out; iz <= _nzhi_out; iz++)
    for (int iy = _nylo_out; iy <= _nyhi_out; iy++)
      for (int ix = _nxlo_out; ix <= _nxhi_out; ix++) {
        numtyp2 v;
        v.x = host_grid_brick[iz][iy][ix];
        v.y = (numtyp)0;
        _cgrid_brick[n] = v;
        n++;
      }
  _cgrid_brick.update_device(_num_grid_points, false);

  _felec = felec;
  fphi_mpole();

  _fdip_sum_phi.update_host(_max_thetai_size*20, false);

  *host_fphi = _fdip_sum_phi.host.begin();
}

// ---------------------------------------------------------------------------
// Interpolate the potential from the PME grid
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int BaseAmoebaT::fphi_mpole() {
  int ainum=ans->inum();
  if (ainum == 0)
    return 0;

  // Compute the block size and grid size to keep all cores busy

  const int BX=block_size();
  const int GX=static_cast<int>(ceil(static_cast<double>(ainum)/BX));

  time_pair.start();
  int ngridxy = _ngridx * _ngridy;
  k_fphi_mpole.set_size(GX,BX);
  k_fphi_mpole.run(&_thetai1, &_thetai2, &_thetai3, &_igrid, &_cgrid_brick,
                  &_fdip_sum_phi, &_bsorder, &ainum, &_felec,
                  &_nzlo_out, &_nylo_out, &_nxlo_out, &ngridxy, &_ngridx);
  time_pair.stop();

  return GX;
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary, and then compute polar real-space
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseAmoebaT::compute_polar_real(int *host_amtype, int *host_amgroup,
                                     double **host_rpole, double **host_uind,
                                     double **host_uinp, double *host_pval,
                                     const bool eflag_in, const bool vflag_in,
                                     const bool eatom, const bool vatom,
                                     const double aewald, const double felec,
                                     const double off2_polar, void **tep_ptr) {

  // cast necessary data arrays from host to device

  cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, host_pval);
  atom->add_extra_data();

  *tep_ptr=_tep.host.begin();

  _off2_polar = off2_polar;
  _felec = felec;
  _aewald = aewald;
  const int red_blocks=polar_real(_eflag,_vflag);

  // only copy answers (forces, energies and virial) back from the device
  //   in the last kernel (which is polar_real here)
  ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  device->add_ans_object(ans);

  // copy tep from device to host
  _tep.update_host(_max_tep_size*4,false);
}

// ---------------------------------------------------------------------------
// Return the memory bytes allocated on the host and device
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
double BaseAmoebaT::host_memory_usage_atomic() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(BaseAmoeba<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Setup the FFT plan: only placeholder for now
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
void BaseAmoebaT::setup_fft(const int /*numel*/, const int /*element_type*/)
{
  // TODO: setting up FFT plan based on the backend (cuFFT or hipFFT)
}

// ---------------------------------------------------------------------------
// Compute FFT on the device: only placeholder for now
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
void BaseAmoebaT::compute_fft1d(void * /*in*/, void * /*out*/,
                                const int /*numel*/, const int /*mode*/)
{
  // TODO: setting up FFT plan based on the backend (cuFFT or hipFFT)
  #if 0 // !defined(USE_OPENCL) && !defined(USE_HIP)
  if (fft_plan_created == false) {
    int m = numel/2;
    cufftPlan1d(&plan, m, CUFFT_Z2Z, 1);
    fft_plan_created = true;
  }

  // n = number of double complex
  int n = numel/2;

  // copy the host array to the device (data)
  UCL_Vector<cufftDoubleComplex,cufftDoubleComplex> data;
  data.alloc(n, *(this->ucl_device), UCL_READ_WRITE, UCL_READ_WRITE);
  int m = 0;
  double* d_in = (double*)in;
  for (int i = 0; i < n; i++) {
    data[i].x = d_in[m];
    data[i].y = d_in[m+1];
    m += 2;
  }
  data.update_device(false);

  // perform the in-place forward FFT

  cufftResult result = cufftExecZ2Z(plan, (cufftDoubleComplex*)&data.device,
    (cufftDoubleComplex*)&data.device, CUFFT_FORWARD);
  if (result != CUFFT_SUCCESS) printf("failed cufft %d\n", result);
  ucl_device->sync();
  data.update_host(false);

  // copy back the data to the host array

  m = 0;
  double* d_out = (double*)out;
  for (int i = 0; i < n; i++) {
    d_out[m] = data[i].x;
    d_out[m+1] = data[i].y;
    m += 2;
  }

  data.clear();
  #endif
}

// ---------------------------------------------------------------------------
// Copy the extra data from host to device
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
void BaseAmoebaT::cast_extra_data(int* amtype, int* amgroup, double** rpole,
                                  double** uind, double** uinp, double* pval) {
  // signal that we need to transfer extra data from the host

  atom->extra_data_unavail();

  int _nall=atom->nall();
  numtyp4 *pextra=reinterpret_cast<numtyp4*>(&(atom->extra[0]));

  int n = 0;
  int nstride = 1; //4;
  if (rpole) {
    for (int i = 0; i < _nall; i++) {
      int idx = n+i*nstride;
      pextra[idx].x = rpole[i][0];
      pextra[idx].y = rpole[i][1];
      pextra[idx].z = rpole[i][2];
      pextra[idx].w = rpole[i][3];
    }

    n += nstride*_nall;
    for (int i = 0; i < _nall; i++) {
      int idx = n+i*nstride;
      pextra[idx].x = rpole[i][4];
      pextra[idx].y = rpole[i][5];
      pextra[idx].z = rpole[i][6];
      pextra[idx].w = rpole[i][8];
    }

    n += nstride*_nall;
    for (int i = 0; i < _nall; i++) {
      int idx = n+i*nstride;
      pextra[idx].x = rpole[i][9];
      pextra[idx].y = rpole[i][12];
      pextra[idx].z = (numtyp)amtype[i];
      pextra[idx].w = (numtyp)amgroup[i];
    }
  } else {
    n += 2*nstride*_nall;
  }

  n += nstride*_nall;
  if (uind) {
    for (int i = 0; i < _nall; i++) {
      int idx = n+i*nstride;
      pextra[idx].x = uind[i][0];
      pextra[idx].y = uind[i][1];
      pextra[idx].z = uind[i][2];
      pextra[idx].w = 0;
    }
  }

  n += nstride*_nall;
  if (uinp) {
    for (int i = 0; i < _nall; i++) {
      int idx = n+i*nstride;
      pextra[idx].x = uinp[i][0];
      pextra[idx].y = uinp[i][1];
      pextra[idx].z = uinp[i][2];
      pextra[idx].w = 0;
    }
  }

  n += nstride*_nall;
  if (pval) {
    for (int i = 0; i < _nall; i++) {
      int idx = n+i*nstride;
      pextra[idx].x = pval[i];
      pextra[idx].y = 0;
      pextra[idx].z = 0;
      pextra[idx].w = 0;
    }
  }
}

// ---------------------------------------------------------------------------
// Compile (load) the kernel strings and set the kernels
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
void BaseAmoebaT::compile_kernels(UCL_Device &dev, const void *pair_str,
                                  const char *kname_multipole,
                                  const char *kname_udirect2b,
                                  const char *kname_umutual2b,
                                  const char *kname_polar,
                                  const char *kname_fphi_uind,
                                  const char *kname_fphi_mpole,
                                  const char *kname_short_nbor,
                                  const char* kname_special15) {
  if (_compiled)
    return;

  if (pair_program) delete pair_program;
  pair_program=new UCL_Program(dev);
  std::string oclstring = device->compile_string()+" -DEVFLAG=1";
  pair_program->load_string(pair_str, oclstring.c_str(),nullptr, screen);

  k_multipole.set_function(*pair_program, kname_multipole);
  k_udirect2b.set_function(*pair_program, kname_udirect2b);
  k_umutual2b.set_function(*pair_program, kname_umutual2b);
  k_polar.set_function(*pair_program, kname_polar);
  k_fphi_uind.set_function(*pair_program, kname_fphi_uind);
  k_fphi_mpole.set_function(*pair_program, kname_fphi_mpole);
  k_short_nbor.set_function(*pair_program, kname_short_nbor);
  k_special15.set_function(*pair_program, kname_special15);
  pos_tex.get_texture(*pair_program, "pos_tex");
  q_tex.get_texture(*pair_program, "q_tex");

  _compiled=true;

  #if defined(USE_OPENCL) && (defined(CL_VERSION_2_1) || defined(CL_VERSION_3_0))
  if (dev.has_subgroup_support()) {
    int mx_subgroup_sz = k_polar.max_subgroup_size(_block_size);
    if (_threads_per_atom > mx_subgroup_sz)
      _threads_per_atom = mx_subgroup_sz;
    device->set_simd_size(mx_subgroup_sz);
  }
  #endif

}

// ---------------------------------------------------------------------------
//  Specify 1-5 neighbors from the current neighbor list
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
int BaseAmoebaT::add_onefive_neighbors() {
  // Compute the block size and grid size to keep all cores busy
  const int BX=block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(ans->inum())/
                               (BX/_threads_per_atom)));

  int _nall=atom->nall();
  int ainum=ans->inum();
  int nbor_pitch=nbor->nbor_pitch();

  k_special15.set_size(GX,BX);
  k_special15.run(&nbor->dev_nbor, &_nbor_data->begin(),
                  &atom->dev_tag, &dev_nspecial15, &dev_special15,
                  &ainum, &_nall, &nbor_pitch,
                  &_threads_per_atom);

  return GX;
}

template class BaseAmoeba<PRECISION,ACC_PRECISION>;
}
