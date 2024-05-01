/***************************************************************************
                                gayberne.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Host code for Gay-Berne potential acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "gayberne_cl.h"
#include "gayberne_lj_cl.h"
#elif defined(USE_CUDART)
const char *gayberne=0;
const char *gayberne_lj=0;
#else
#include "gayberne_cubin.h"
#include "gayberne_lj_cubin.h"
#endif

#include "lal_gayberne.h"
#include <cassert>
namespace LAMMPS_AL {

#define GayBerneT GayBerne<numtyp, acctyp>
extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
GayBerneT::GayBerne() : BaseEllipsoid<numtyp,acctyp>(),
                                  _allocated(false) {
}

template <class numtyp, class acctyp>
GayBerneT::~GayBerne() {
  clear();
}

template <class numtyp, class acctyp>
int GayBerneT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_ellipsoid(max_nbors);
}

template <class numtyp, class acctyp>
int GayBerneT::init(const int ntypes, const double gamma,
                         const double upsilon, const double mu,
                         double **host_shape, double **host_well,
                         double **host_cutsq, double **host_sigma,
                         double **host_epsilon, double *host_lshape,
                         int **h_form, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4,
                         double **host_offset, const double *host_special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size,
                         const double gpu_split, FILE *_screen) {
  int success;
  success=this->init_base(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                          _screen,ntypes,h_form,gayberne,gayberne_lj,
                          "k_gayberne");
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  _shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->block_size()>=max_shared_types) {
    lj_types=max_shared_types;
    _shared_types=true;
  }
  _lj_types=lj_types;

  // Allocate a host write buffer for copying type data
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write[i]=0.0;

  sigma_epsilon.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,lj_types,sigma_epsilon,host_write,
                         host_sigma,host_epsilon);

  this->cut_form.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,lj_types,this->cut_form,host_write,
                         host_cutsq,h_form);

  lj1.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj1,host_write,host_lj1,host_lj2,
                         host_cutsq,h_form);

  lj3.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj3,host_write,host_lj3,host_lj4,
                         host_offset);

  dev_error.alloc(1,*(this->ucl_device),UCL_WRITE_ONLY);
  dev_error.zero();

  // Allocate, cast and asynchronous memcpy of constant data
  // Copy data for bonded interactions
  gamma_upsilon_mu.alloc(7,*(this->ucl_device),UCL_READ_ONLY);
  host_write[0]=static_cast<numtyp>(gamma);
  host_write[1]=static_cast<numtyp>(upsilon);
  host_write[2]=static_cast<numtyp>(mu);
  host_write[3]=static_cast<numtyp>(host_special_lj[0]);
  host_write[4]=static_cast<numtyp>(host_special_lj[1]);
  host_write[5]=static_cast<numtyp>(host_special_lj[2]);
  host_write[6]=static_cast<numtyp>(host_special_lj[3]);
  ucl_copy(gamma_upsilon_mu,host_write,7,false);

  lshape.alloc(ntypes,*(this->ucl_device),UCL_READ_ONLY);
  UCL_H_Vec<double> d_view;
  d_view.view(host_lshape,lshape.numel(),*(this->ucl_device));
  ucl_copy(lshape,d_view,false);

  // Copy shape, well, sigma, epsilon, and cutsq onto GPU
  // - cast if necessary
  shape.alloc(ntypes,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<ntypes; i++) {
    host_write[i*4]=host_shape[i][0];
    host_write[i*4+1]=host_shape[i][1];
    host_write[i*4+2]=host_shape[i][2];
  }
  UCL_H_Vec<numtyp4> view4;
  view4.view(host_write,shape.numel());
  ucl_copy(shape,view4,false);

  well.alloc(ntypes,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<ntypes; i++) {
    host_write[i*4]=host_well[i][0];
    host_write[i*4+1]=host_well[i][1];
    host_write[i*4+2]=host_well[i][2];
  }
  view4.view(host_write,well.numel());
  ucl_copy(well,view4,false);

  _allocated=true;
  this->_max_bytes=sigma_epsilon.row_bytes()+this->cut_form.row_bytes()+
                   lj1.row_bytes()+lj3.row_bytes()+gamma_upsilon_mu.row_bytes()+
                   lshape.row_bytes()+shape.row_bytes()+well.row_bytes();

  return 0;
}

template <class numtyp, class acctyp>
void GayBerneT::clear() {
  if (!_allocated)
    return;

  UCL_H_Vec<int> err_flag(1,*(this->ucl_device));
  ucl_copy(err_flag,dev_error,false);
  if (err_flag[0] == 2)
    std::cerr << "BAD MATRIX INVERSION IN FORCE COMPUTATION.\n";
  err_flag.clear();

  _allocated=false;

  dev_error.clear();
  lj1.clear();
  lj3.clear();
  sigma_epsilon.clear();
  this->cut_form.clear();

  shape.clear();
  well.clear();
  lshape.clear();
  gamma_upsilon_mu.clear();

  this->clear_base();
}

template <class numtyp, class acctyp>
double GayBerneT::host_memory_usage() const {
  return this->host_memory_usage_base()+sizeof(GayBerneT)+
         4*sizeof(numtyp);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int GayBerneT::loop(const int eflag, const int vflag) {
  const int BX=this->block_size();
  int GX=0, NGX;
  int stride=this->nbor->nbor_pitch();
  int ainum=this->ans->inum();

  if (this->_multiple_forms) {
    this->time_nbor1.start();
    if (this->_last_ellipse>0) {
      // ------------ ELLIPSE_ELLIPSE and ELLIPSE_SPHERE ---------------
      GX=static_cast<int>(ceil(static_cast<double>(this->_last_ellipse)/
                               (BX/this->_threads_per_atom)));
      NGX=static_cast<int>(ceil(static_cast<double>(this->_last_ellipse)/BX));
      this->pack_nbors(NGX,BX, 0, this->_last_ellipse,ELLIPSE_SPHERE,
                                         ELLIPSE_ELLIPSE,_shared_types,_lj_types);
      this->time_nbor1.stop();

      this->time_ellipsoid.start();
      this->k_elps_sel->set_size(GX,BX);
      this->k_elps_sel->run(&this->atom->x, &this->atom->quat,
                            &this->shape, &this->well, &this->gamma_upsilon_mu,
                            &this->sigma_epsilon, &this->_lj_types,
                            &this->lshape, &this->nbor->dev_nbor, &stride,
                            &this->ans->force, &ainum, &this->ans->engv,
                            &this->dev_error, &eflag, &vflag,
                            &this->_last_ellipse, &this->_threads_per_atom);
      this->time_ellipsoid.stop();

      if (this->_last_ellipse==this->ans->inum()) {
        this->time_nbor2.start();
        this->time_nbor2.stop();
        this->time_ellipsoid2.start();
        this->time_ellipsoid2.stop();
        this->time_lj.start();
        this->time_lj.stop();
        return ainum;
      }

      // ------------ SPHERE_ELLIPSE ---------------

      this->time_nbor2.start();
      GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum()-
                               this->_last_ellipse)/
                               (BX/this->_threads_per_atom)));
      NGX=static_cast<int>(ceil(static_cast<double>(this->ans->inum()-
                                this->_last_ellipse)/BX));
      this->pack_nbors(NGX,BX,this->_last_ellipse,this->ans->inum(),
                                         SPHERE_ELLIPSE,SPHERE_ELLIPSE,_shared_types,_lj_types);
      this->time_nbor2.stop();

      this->time_ellipsoid2.start();
      this->k_sphere_elps_sel->set_size(GX,BX);
      this->k_sphere_elps_sel->run(&this->atom->x, &this->atom->quat,
                                   &this->shape,  &this->well,
                                   &this->gamma_upsilon_mu,
                                   &this->sigma_epsilon, &this->_lj_types,
                                   &this->lshape,  &this->nbor->dev_nbor,
                                   &stride, &this->ans->force,
                                   &this->ans->engv, &this->dev_error,
                                   &eflag, &vflag, &this->_last_ellipse,
                                   &ainum, &this->_threads_per_atom);
      this->time_ellipsoid2.stop();
   } else {
      GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum()-
                               this->_last_ellipse)/
                               (BX/this->_threads_per_atom)));
      this->ans->force.zero();
      this->ans->engv.zero();
      this->time_nbor1.stop();
      this->time_ellipsoid.start();
      this->time_ellipsoid.stop();
      this->time_nbor2.start();
      this->time_nbor2.stop();
      this->time_ellipsoid2.start();
      this->time_ellipsoid2.stop();
    }

    // ------------         LJ      ---------------
    this->time_lj.start();
    if (this->_last_ellipse<this->ans->inum()) {
      if (this->_shared_types) {
        this->k_lj_sel->set_size(GX,BX);
        this->k_lj_sel->run(&this->atom->x, &this->lj1, &this->lj3,
                            &this->gamma_upsilon_mu, &stride,
                            &this->nbor->dev_packed, &this->ans->force,
                            &this->ans->engv, &this->dev_error, &eflag,
                            &vflag, &this->_last_ellipse, &ainum,
                            &this->_threads_per_atom);
      } else {
        this->k_lj.set_size(GX,BX);
        this->k_lj.run(&this->atom->x, &this->lj1, &this->lj3,
                       &this->_lj_types, &this->gamma_upsilon_mu, &stride,
                       &this->nbor->dev_packed, &this->ans->force,
                       &this->ans->engv, &this->dev_error, &eflag,
                       &vflag, &this->_last_ellipse, &ainum,
                       &this->_threads_per_atom);
      }
    }
    this->time_lj.stop();
  } else {
    GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                             (BX/this->_threads_per_atom)));
    NGX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/BX));
    this->time_nbor1.start();
    this->pack_nbors(NGX, BX, 0, this->ans->inum(),SPHERE_SPHERE,
                                 ELLIPSE_ELLIPSE,_shared_types,_lj_types);
    this->time_nbor1.stop();
    this->time_ellipsoid.start();
    this->k_elps_sel->set_size(GX,BX);
    this->k_elps_sel->run(&this->atom->x,  &this->atom->quat,
                          &this->shape, &this->well, &this->gamma_upsilon_mu,
                          &this->sigma_epsilon, &this->_lj_types, &this->lshape,
                          &this->nbor->dev_nbor, &stride, &this->ans->force,
                          &ainum,  &this->ans->engv, &this->dev_error,
                          &eflag, &vflag, &ainum, &this->_threads_per_atom);
    this->time_ellipsoid.stop();
  }
  return ainum;
}

template class GayBerne<PRECISION,ACC_PRECISION>;
}
