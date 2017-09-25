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
#include "error.h"
#include "force.h"
#include "kspace.h"

using namespace LAMMPS_NS;

#define BIG   1.0e20
#define SMALL 1.0e-4

/* ---------------------------------------------------------------------- */

DomainKokkos::DomainKokkos(LAMMPS *lmp) : Domain(lmp) {}

/* ---------------------------------------------------------------------- */

void DomainKokkos::init()
{
  atomKK = (AtomKokkos *) atom;
  Domain::init();
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxlo/hi
   for triclinic, atoms must be in lamda coords (0-1) before reset_box is called
------------------------------------------------------------------------- */

template<class DeviceType>
struct DomainResetBoxFunctor{
public:
  typedef DeviceType device_type;
  typename ArrayTypes<DeviceType>::t_x_array x;

  struct value_type {
    double value[3][2] ;
  };

  DomainResetBoxFunctor(DAT::tdual_x_array _x):
    x(_x.view<DeviceType>()) {}

  KOKKOS_INLINE_FUNCTION
  void init(value_type &dst) const {
    dst.value[2][0] = dst.value[1][0] = dst.value[0][0] = BIG;
    dst.value[2][1] = dst.value[1][1] = dst.value[0][1] = -BIG;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type &dst,
             const volatile value_type &src) const {
    dst.value[0][0] = MIN(dst.value[0][0],src.value[0][0]);
    dst.value[0][1] = MAX(dst.value[0][1],src.value[0][1]);
    dst.value[1][0] = MIN(dst.value[1][0],src.value[1][0]);
    dst.value[1][1] = MAX(dst.value[1][1],src.value[1][1]);
    dst.value[2][0] = MIN(dst.value[2][0],src.value[2][0]);
    dst.value[2][1] = MAX(dst.value[2][1],src.value[2][1]);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i, value_type &dst) const {
    dst.value[0][0] = MIN(dst.value[0][0],x(i,0));
    dst.value[0][1] = MAX(dst.value[0][1],x(i,0));
    dst.value[1][0] = MIN(dst.value[1][0],x(i,1));
    dst.value[1][1] = MAX(dst.value[1][1],x(i,1));
    dst.value[2][0] = MIN(dst.value[2][0],x(i,2));
    dst.value[2][1] = MAX(dst.value[2][1],x(i,2));
  }
};

void DomainKokkos::reset_box()
{
  // perform shrink-wrapping
  // compute extent of atoms on this proc
  // for triclinic, this is done in lamda space

  atomKK->sync(Device,X_MASK);

  if (nonperiodic == 2) {

    int nlocal = atom->nlocal;

    DomainResetBoxFunctor<LMPDeviceType>::value_type result;

    DomainResetBoxFunctor<LMPDeviceType>
      f(atomKK->k_x);
    Kokkos::parallel_reduce(nlocal,f,result);

    double (*extent)[2] = result.value;
    double all[3][2];

    // compute extent across all procs
    // flip sign of MIN to do it in one Allreduce MAX

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);

    // for triclinic, convert back to box coords before changing box

    if (triclinic) lamda2x(atom->nlocal);

    // in shrink-wrapped dims, set box by atom extent
    // if minimum set, enforce min box size settings
    // for triclinic, convert lamda extent to box coords, then set box lo/hi
    // decided NOT to do the next comment - don't want to sneakily change tilt
    // for triclinic, adjust tilt factors if 2nd dim is shrink-wrapped,
    //   so that displacement in 1st dim stays the same

    if (triclinic == 0) {
      if (xperiodic == 0) {
        if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - small[0];
        else if (boundary[0][0] == 3)
          boxlo[0] = MIN(-all[0][0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = all[0][1] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(all[0][1]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - small[1];
        else if (boundary[1][0] == 3)
          boxlo[1] = MIN(-all[1][0]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = all[1][1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(all[1][1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
      }
      if (zperiodic == 0) {
        if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - small[2];
        else if (boundary[2][0] == 3)
          boxlo[2] = MIN(-all[2][0]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = all[2][1] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(all[2][1]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
      }

    } else {
      double lo[3],hi[3];
      if (xperiodic == 0) {
        lo[0] = -all[0][0]; lo[1] = 0.0; lo[2] = 0.0;
        Domain::lamda2x(lo,lo);
        hi[0] = all[0][1]; hi[1] = 0.0; hi[2] = 0.0;
        Domain::lamda2x(hi,hi);
        if (boundary[0][0] == 2) boxlo[0] = lo[0] - small[0];
        else if (boundary[0][0] == 3) boxlo[0] = MIN(lo[0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = hi[0] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(hi[0]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        lo[0] = 0.0; lo[1] = -all[1][0]; lo[2] = 0.0;
        Domain::lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = all[1][1]; hi[2] = 0.0;
        Domain::lamda2x(hi,hi);
        if (boundary[1][0] == 2) boxlo[1] = lo[1] - small[1];
        else if (boundary[1][0] == 3) boxlo[1] = MIN(lo[1]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = hi[1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(hi[1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
        //xy *= (boxhi[1]-boxlo[1]) / yprd;
      }
      if (zperiodic == 0) {
        lo[0] = 0.0; lo[1] = 0.0; lo[2] = -all[2][0];
        Domain::lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = 0.0; hi[2] = all[2][1];
        Domain::lamda2x(hi,hi);
        if (boundary[2][0] == 2) boxlo[2] = lo[2] - small[2];
        else if (boundary[2][0] == 3) boxlo[2] = MIN(lo[2]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = hi[2] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(hi[2]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
        //xz *= (boxhi[2]-boxlo[2]) / xprd;
        //yz *= (boxhi[2]-boxlo[2]) / yprd;
      }
    }
  }

  // reset box whether shrink-wrapping or not

  set_global_box();
  set_local_box();

  // if shrink-wrapped & kspace is defined (i.e. using MSM), call setup()
  // also call init() (to test for compatibility) ?

  if (nonperiodic == 2 && force->kspace) {
    //force->kspace->init();
    force->kspace->setup();
  }

  // if shrink-wrapped & triclinic, re-convert to lamda coords for new box
  // re-invoke pbc() b/c x2lamda result can be outside [0,1] due to roundoff

  if (nonperiodic == 2 && triclinic) {
    x2lamda(atom->nlocal);
    pbc();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, int PERIODIC, int DEFORM_VREMAP>
struct DomainPBCFunctor {
  typedef DeviceType device_type;
  double lo[3],hi[3],period[3];
  typename ArrayTypes<DeviceType>::t_x_array x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_int_1d mask;
  typename ArrayTypes<DeviceType>::t_imageint_1d image;
  int deform_groupbit;
  double h_rate[6];
  int xperiodic,yperiodic,zperiodic;

  DomainPBCFunctor(double* _lo, double* _hi, double* _period,
                   DAT::tdual_x_array _x, DAT::tdual_v_array _v,
                   DAT::tdual_int_1d _mask, DAT::tdual_imageint_1d _image,
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
        imageint idim = image[i] & IMGMASK;
        const imageint otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x(i,0) >= hi[0]) {
        x(i,0) -= period[0];
        x(i,0) = MAX(x(i,0),lo[0]);
        if (DEFORM_VREMAP && (mask[i] & deform_groupbit)) v(i,0) -= h_rate[0];
        imageint idim = image[i] & IMGMASK;
        const imageint otherdims = image[i] ^ idim;
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
        imageint idim = (image[i] >> IMGBITS) & IMGMASK;
        const imageint otherdims = image[i] ^ (idim << IMGBITS);
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
        imageint idim = (image[i] >> IMGBITS) & IMGMASK;
        const imageint otherdims = image[i] ^ (idim << IMGBITS);
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
        imageint idim = image[i] >> IMG2BITS;
        const imageint otherdims = image[i] ^ (idim << IMG2BITS);
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
        imageint idim = image[i] >> IMG2BITS;
        const imageint otherdims = image[i] ^ (idim << IMG2BITS);
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

  atomKK->modified(Device,X_MASK|V_MASK|IMAGE_MASK);
}

/* ----------------------------------------------------------------------
   remap all points into the periodic box no matter how far away
   adjust 3 image flags encoded in image accordingly
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before doing remap
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void DomainKokkos::remap_all()
{
  atomKK->sync(Device,X_MASK | IMAGE_MASK);

  x = atomKK->k_x.view<LMPDeviceType>();
  image = atomKK->k_image.view<LMPDeviceType>();
  int nlocal = atomKK->nlocal;

  if (triclinic == 0) {
    for (int i=0; i<3; i++) {
      lo[i] = boxlo[i];
      hi[i] = boxhi[i];
      period[i] = prd[i];
    }
  } else {
    for (int i=0; i<3; i++) {
      lo[i] = boxlo_lamda[i];
      hi[i] = boxhi_lamda[i];
      period[i] = prd_lamda[i];
    }
    x2lamda(nlocal);
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagDomain_remap_all>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(Device,X_MASK | IMAGE_MASK);

  if (triclinic) lamda2x(nlocal);
}

KOKKOS_INLINE_FUNCTION
void DomainKokkos::operator()(TagDomain_remap_all, const int &i) const {
    imageint idim,otherdims;
    if (xperiodic) {
      while (x(i,0) < lo[0]) {
        x(i,0) += period[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      while (x(i,0) >= hi[0]) {
        x(i,0) -= period[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      x(i,0) = MAX(x(i,0),lo[0]);
    }

    if (yperiodic) {
      while (x(i,1) < lo[1]) {
        x(i,1) += period[1];
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      while (x(i,1) >= hi[1]) {
        x(i,1) -= period[1];
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      x(i,1) = MAX(x(i,1),lo[1]);
    }

    if (zperiodic) {
      while (x(i,2) < lo[2]) {
        x(i,2) += period[2];
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      while (x(i,2) >= hi[2]) {
        x(i,2) -= period[2];
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      x(i,2) = MAX(x(i,2),lo[2]);
    }
}

/* ----------------------------------------------------------------------
   adjust image flags due to triclinic box flip
   flip operation is changing box vectors A,B,C to new A',B',C'
     A' = A              (A does not change)
     B' = B + mA         (B shifted by A)
     C' = C + pB + nA    (C shifted by B and/or A)
   this requires the image flags change from (a,b,c) to (a',b',c')
   so that x_unwrap for each atom is same before/after
     x_unwrap_before = xlocal + aA + bB + cC
     x_unwrap_after = xlocal + a'A' + b'B' + c'C'
   this requires:
     c' = c
     b' = b - cp
     a' = a - (b-cp)m - cn = a - b'm - cn
   in other words, for xy flip, change in x flag depends on current y flag
   this is b/c the xy flip dramatically changes which tiled image of
     simulation box an unwrapped point maps to
------------------------------------------------------------------------- */

void DomainKokkos::image_flip(int m_in, int n_in, int p_in)
{
  m_flip = m_in;
  n_flip = n_in;
  p_flip = p_in;

  atomKK->sync(Device,IMAGE_MASK);

  image = atomKK->k_image.view<LMPDeviceType>();
  int nlocal = atomKK->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagDomain_image_flip>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(Device,IMAGE_MASK);
}

KOKKOS_INLINE_FUNCTION
void DomainKokkos::operator()(TagDomain_image_flip, const int &i) const {
  int xbox = (image[i] & IMGMASK) - IMGMAX;
  int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image[i] >> IMG2BITS) - IMGMAX;

  ybox -= p_flip*zbox;
  xbox -= m_flip*ybox + n_flip*zbox;

  image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
    (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
    (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for all N atoms
   x = H lamda + x0;
------------------------------------------------------------------------- */

void DomainKokkos::lamda2x(int n)
{
  atomKK->sync(Device,X_MASK);

  x = atomKK->k_x.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagDomain_lamda2x>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,X_MASK);
}

KOKKOS_INLINE_FUNCTION
void DomainKokkos::operator()(TagDomain_lamda2x, const int &i) const {
  x(i,0) = h[0]*x(i,0) + h[5]*x(i,1) + h[4]*x(i,2) + boxlo[0];
  x(i,1) = h[1]*x(i,1) + h[3]*x(i,2) + boxlo[1];
  x(i,2) = h[2]*x(i,2) + boxlo[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for all N atoms
   lamda = H^-1 (x - x0)
------------------------------------------------------------------------- */

void DomainKokkos::x2lamda(int n)
{
  atomKK->sync(Device,X_MASK);

  x = atomKK->k_x.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagDomain_x2lamda>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,X_MASK);
}

KOKKOS_INLINE_FUNCTION
void DomainKokkos::operator()(TagDomain_x2lamda, const int &i) const {
  F_FLOAT delta[3];
  delta[0] = x(i,0) - boxlo[0];
  delta[1] = x(i,1) - boxlo[1];
  delta[2] = x(i,2) - boxlo[2];

  x(i,0) = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  x(i,1) = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  x(i,2) = h_inv[2]*delta[2];
}

