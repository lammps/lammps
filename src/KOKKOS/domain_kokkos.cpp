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

#include "domain_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DomainKokkos::DomainKokkos(LAMMPS *lmp) : Domain(lmp) {}

/* ---------------------------------------------------------------------- */

void DomainKokkos::init()
{
  atomKK = (AtomKokkos *) atom;
  Domain::init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, int PERIODIC, int DEFORM_VREMAP>
struct DomainPBCFunctor {
  typedef DeviceType device_type;
  double lo[3],hi[3],period[3];
  typename ArrayTypes<DeviceType>::t_x_array x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_int_1d mask;
  typename ArrayTypes<DeviceType>::t_int_1d image;
  int deform_groupbit;
  double h_rate[6];
  int xperiodic,yperiodic,zperiodic;

  DomainPBCFunctor(double* _lo, double* _hi, double* _period,
                   DAT::tdual_x_array _x, DAT::tdual_v_array _v,
                   DAT::tdual_int_1d _mask, DAT::tdual_int_1d _image, 
                   int _deform_groupbit, double* _h_rate,
                   int _xperiodic, int _yperiodic, int _zperiodic):
    x(_x.view<DeviceType>()), v(_v.view<DeviceType>()),
    mask(_mask.view<DeviceType>()), image(_image.view<DeviceType>()),
    deform_groupbit(_deform_groupbit),
    xperiodic(_xperiodic), yperiodic(_yperiodic), zperiodic(_zperiodic){
    lo[0]=_lo[0]; lo[1]=_lo[1]; lo[2]=_lo[2];
    hi[0]=_hi[0]; hi[1]=_hi[1]; hi[2]=_hi[2];
    period[0]=_period[0]; period[1]=_period[1]; period[2]=_period[2];
    h_rate[0]=_h_rate[0]; h_rate[1]=_h_rate[1]; h_rate[2]=_h_rate[2];
    h_rate[3]=_h_rate[3]; h_rate[4]=_h_rate[4]; h_rate[5]=_h_rate[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const {
    if (PERIODIC && xperiodic) {
      if (x(i,0) < lo[0]) {
        x(i,0) += period[0];
        if (DEFORM_VREMAP && (mask[i] & deform_groupbit)) v(i,0) += h_rate[0];
        int idim = image[i] & IMGMASK;
        const int otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x(i,0) >= hi[0]) {
        x(i,0) -= period[0];
        x(i,0) = MAX(x(i,0),lo[0]);
        if (DEFORM_VREMAP && (mask[i] & deform_groupbit)) v(i,0) -= h_rate[0];
        int idim = image[i] & IMGMASK;
        const int otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
    }
    
    if (PERIODIC && yperiodic) {
      if (x(i,1) < lo[1]) {
        x(i,1) += period[1];
        if (DEFORM_VREMAP && (mask[i] & deform_groupbit)) {
          v(i,0) += h_rate[5];
          v(i,1) += h_rate[1];
        }
        int idim = (image[i] >> IMGBITS) & IMGMASK;
        const int otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      if (x(i,1) >= hi[1]) {
        x(i,1) -= period[1];
        x(i,1) = MAX(x(i,1),lo[1]);
        if (DEFORM_VREMAP && (mask[i] & deform_groupbit)) {
          v(i,0) -= h_rate[5];
          v(i,1) -= h_rate[1];
        }
        int idim = (image[i] >> IMGBITS) & IMGMASK;
        const int otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
    }
    
    if (PERIODIC && zperiodic) {
      if (x(i,2) < lo[2]) {
        x(i,2) += period[2];
        if (DEFORM_VREMAP && (mask[i] & deform_groupbit)) {
          v(i,0) += h_rate[4];
          v(i,1) += h_rate[3];
          v(i,2) += h_rate[2];
        }
        int idim = image[i] >> IMG2BITS;
        const int otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      if (x(i,2) >= hi[2]) {
        x(i,2) -= period[2];
        x(i,2) = MAX(x(i,2),lo[2]);
        if (DEFORM_VREMAP && (mask[i] & deform_groupbit)) {
          v(i,0) -= h_rate[4];
          v(i,1) -= h_rate[3];
          v(i,2) -= h_rate[2];
        }
        int idim = image[i] >> IMG2BITS;
        const int otherdims = image[i] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
    }
  }
};

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   if fix deform, remap velocity of fix group atoms by box edge velocities
   for triclinic, atoms must be in lamda coords (0-1) before pbc is called
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void DomainKokkos::pbc()
{
  double *lo,*hi,*period;
  int nlocal = atomKK->nlocal;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }

  atomKK->sync(Device,X_MASK|V_MASK|MASK_MASK|IMAGE_MASK);
  atomKK->modified(Device,X_MASK|V_MASK);

  if (xperiodic || yperiodic || zperiodic) {
    if (deform_vremap) {
      DomainPBCFunctor<LMPDeviceType,1,1> 
        f(lo,hi,period,
          atomKK->k_x,atomKK->k_v,atomKK->k_mask,atomKK->k_image,
          deform_groupbit,h_rate,xperiodic,yperiodic,zperiodic);
      Kokkos::parallel_for(nlocal,f);
    } else {
      DomainPBCFunctor<LMPDeviceType,1,0> 
        f(lo,hi,period,
          atomKK->k_x,atomKK->k_v,atomKK->k_mask,atomKK->k_image,
          deform_groupbit,h_rate,xperiodic,yperiodic,zperiodic);
      Kokkos::parallel_for(nlocal,f);
    }
  } else {
    if (deform_vremap) {
      DomainPBCFunctor<LMPDeviceType,0,1> 
        f(lo,hi,period,
          atomKK->k_x,atomKK->k_v,atomKK->k_mask,atomKK->k_image,
          deform_groupbit,h_rate,xperiodic,yperiodic,zperiodic);
      Kokkos::parallel_for(nlocal,f);
    } else {
      DomainPBCFunctor<LMPDeviceType,0,0> 
        f(lo,hi,period,
          atomKK->k_x,atomKK->k_v,atomKK->k_mask,atomKK->k_image,
          deform_groupbit,h_rate,xperiodic,yperiodic,zperiodic);
      Kokkos::parallel_for(nlocal,f);
    }
  }

  LMPDeviceType::fence();
}

