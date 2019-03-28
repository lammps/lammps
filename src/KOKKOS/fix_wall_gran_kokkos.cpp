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

#include "fix_wall_gran_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "atom_vec_kokkos.h"
#include "atom_masks.h"
#include "update.h"

using namespace LAMMPS_NS;

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER,REGION};
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY,BONDED_HISTORY};
enum{NONE,CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallGranKokkos<DeviceType>::FixWallGranKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | OMEGA_MASK | TORQUE_MASK | RADIUS_MASK | RMASS_MASK | MASK_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;

  memory->destroy(shearone);
  shearone = NULL;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallGranKokkos<DeviceType>::~FixWallGranKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_shearone, shearone);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::init()
{
  FixWallGran::init();

  if (fix_rigid)
    error->all(FLERR, "wall/gran/kk not yet compatible with rigid.");  
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // do not update shear history during setup

  shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly

  wlo = lo;
  whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      wlo = lo + amplitude - amplitude*cos(arg);
      whi = hi + amplitude - amplitude*cos(arg);
    }
    vwall[axis] = amplitude*omega*sin(arg);
  } else if (wshear) vwall[axis] = vshear;

  copymode = 1;

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega_ = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius_ = atomKK->k_radius.view<DeviceType>();
  int nlocal = atom->nlocal;

  if (pairstyle == HOOKE)
    error->all(FLERR, "wall/gran/kk doesn't yet support hooke style.");
  else if (pairstyle == HOOKE_HISTORY) {
    if (wallstyle == XPLANE) {
      FixWallGranKokkosHookeHistoryFunctor<DeviceType, XPLANE> f(this);
      Kokkos::parallel_for(nlocal,f);
    } else if (wallstyle == YPLANE) {
      FixWallGranKokkosHookeHistoryFunctor<DeviceType, YPLANE> f(this);
      Kokkos::parallel_for(nlocal,f);
    } else if (wallstyle == ZPLANE) {
      FixWallGranKokkosHookeHistoryFunctor<DeviceType, ZPLANE> f(this);
      Kokkos::parallel_for(nlocal,f);
    } else if (wallstyle == ZCYLINDER) {
      FixWallGranKokkosHookeHistoryFunctor<DeviceType, ZCYLINDER> f(this);
      Kokkos::parallel_for(nlocal,f);
    }
  }
  else if (pairstyle == HERTZ_HISTORY)
    error->all(FLERR, "wall/gran/kk doesn't yet support hertz/history style.");

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int WallStyle>
void FixWallGranKokkos<DeviceType>::hooke_history_item(const int &i) const
{
  double vwall_[3];
  vwall_[0] = vwall[0];
  vwall_[1] = vwall[1];
  vwall_[2] = vwall[2];
  
  if (mask[i] & groupbit) {
    X_FLOAT radius = radius_(i);

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    
    if (WallStyle == XPLANE) {
      X_FLOAT del1 = x(i,0) - wlo;
      double del2 = whi - x(i,0);
      if (del1 < del2) dx = del1;
      else dx = -del2;
    } else if (WallStyle == YPLANE) {
      double del1 = x(i,1) - wlo;
      double del2 = whi - x(i,1);
      if (del1 < del2) dy = del1;
      else dy = -del2;
    } else if (WallStyle == ZPLANE) {
      double del1 = x(i,2) - wlo;
      double del2 = whi - x(i,2);
      if (del1 < del2) dz = del1;
      else dz = -del2;
    } else if (WallStyle == ZCYLINDER) {
      double delxy = sqrt(x(i,0)*x(i,0) + x(i,1)*x(i,1));
      double delr = cylradius - delxy;
      if (delr > radius) {
    	dz = cylradius;
      } else {
    	dx = -delr/delxy * x(i,0);
    	dy = -delr/delxy * x(i,1);
     	if (wshear && axis != 2) {
    	  vwall_[0] += vshear * x(i,1)/delxy;
    	  vwall_[1] += -vshear * x(i,0)/delxy;
    	  vwall_[2] = 0.0;
    	}
      }
    }

    double rsq = dx*dx + dy*dy + dz*dz;

    if (rsq > radius*radius) {
      if (history)
    	for (int j = 0; j < 3; j++)
    	  d_shearone(i,j) = 0.0;
    } else {
      // meff = effective mass of sphere
      double meff = rmass(i);
      double r = sqrt(rsq);
      double rinv = 1.0/r;
      double rsqinv = 1.0/rsq;

      // relative translational velocity

      double vr1 = v(i,0) - vwall_[0];
      double vr2 = v(i,1) - vwall_[1];
      double vr3 = v(i,2) - vwall_[2];

      // normal component

      double vnnr = vr1*dx + vr2*dy + vr3*dz;
      double vn1 = dx*vnnr * rsqinv;
      double vn2 = dy*vnnr * rsqinv;
      double vn3 = dz*vnnr * rsqinv;

      // tangential component

      double vt1 = vr1 - vn1;
      double vt2 = vr2 - vn2;
      double vt3 = vr3 - vn3;

      // relative rotational velocity

      double wr1 = radius*omega_(i,0) * rinv;
      double wr2 = radius*omega_(i,1) * rinv;
      double wr3 = radius*omega_(i,2) * rinv;

      // normal forces = Hookian contact + normal velocity damping

      double damp = meff*gamman*vnnr*rsqinv;
      double ccel = kn*(radius-r)*rinv - damp;

      // relative velocities

      double vtr1 = vt1 - (dz*wr2-dy*wr3);
      double vtr2 = vt2 - (dx*wr3-dz*wr1);
      double vtr3 = vt3 - (dy*wr1-dx*wr2);
      double vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
      vrel = sqrt(vrel);

      // shear history effects

      if (shearupdate) {
    	d_shearone(i,0) += vtr1*dt;
    	d_shearone(i,1) += vtr2*dt;
    	d_shearone(i,2) += vtr3*dt;
      }
      double shrmag = sqrt(d_shearone(i,0)*d_shearone(i,0) + d_shearone(i,1)*d_shearone(i,1) + d_shearone(i,2)*d_shearone(i,2));

      // rotate shear displacements

      double rsht = d_shearone(i,0)*dx + d_shearone(i,1)*dy + d_shearone(i,2)*dz;
      rsht = rsht*rsqinv;
      if (shearupdate) {
    	d_shearone(i,0) -= rsht*dx;
    	d_shearone(i,1) -= rsht*dy;
    	d_shearone(i,2) -= rsht*dz;
      }

      // tangential forces = shear + tangential velocity damping

      double fs1 = - (kt*d_shearone(i,0) + meff*gammat*vtr1);
      double fs2 = - (kt*d_shearone(i,1) + meff*gammat*vtr2);
      double fs3 = - (kt*d_shearone(i,2) + meff*gammat*vtr3);

      // rescale frictional displacements and forces if needed

      double fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
      double fn = xmu * fabs(ccel*r);

      if (fs > fn) {
    	if (shrmag != 0.0) {
    	  d_shearone(i,0) = (fn/fs) * (d_shearone(i,0) + meff*gammat*vtr1/kt) -
    	    meff*gammat*vtr1/kt;
    	  d_shearone(i,1) = (fn/fs) * (d_shearone(i,1) + meff*gammat*vtr2/kt) -
    	    meff*gammat*vtr2/kt;
    	  d_shearone(i,2) = (fn/fs) * (d_shearone(i,2) + meff*gammat*vtr3/kt) -
    	    meff*gammat*vtr3/kt;
    	  fs1 *= fn/fs ;
    	  fs2 *= fn/fs;
    	  fs3 *= fn/fs;
    	} else fs1 = fs2 = fs3 = 0.0;
      }

      // forces & torques

      double fx = dx*ccel + fs1;
      double fy = dy*ccel + fs2;
      double fz = dz*ccel + fs3;
      f(i,0) += fx;
      f(i,1) += fy;
      f(i,2) += fz;

      double tor1 = rinv * (dy*fs3 - dz*fs2);
      double tor2 = rinv * (dz*fs1 - dx*fs3);
      double tor3 = rinv * (dx*fs2 - dy*fs1);
      torque(i,0) -= radius*tor1;
      torque(i,1) -= radius*tor2;
      torque(i,2) -= radius*tor3;
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::grow_arrays(int nmax)
{
  if (history) {
    k_shearone.template sync<LMPHostType>(); // force reallocation on host 
    memoryKK->grow_kokkos(k_shearone,shearone,nmax,sheardim,"wall/gran/kk:shearone");
    d_shearone = k_shearone.template view<DeviceType>();
    k_shearone.template modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  if (history) {
    k_shearone.template sync<LMPHostType>();
    for (int m = 0; m < sheardim; m++)
      shearone[j][m] = shearone[i][m];
    k_shearone.template modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int FixWallGranKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_shearone.template sync<LMPHostType>();

  int n = 0;
  for (int j = 0; j < sheardim; j++)
    buf[n++] = shearone[i][j];
  return n;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int FixWallGranKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  for (int j = 0; j < sheardim; j++)
    shearone[nlocal][j] = buf[n++];

  k_shearone.template modify<LMPHostType>();

  return n;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct FixWallGranKokkos_PackExchangeFunctor
{
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  typename AT::t_float_2d _shearone;
  typename AT::t_xfloat_1d_um _buf;
  const int _dnum;

  FixWallGranKokkos_PackExchangeFunctor(
    const typename AT::tdual_xfloat_2d &buf,
    const typename AT::tdual_int_1d &sendlist,
    const typename AT::tdual_int_1d &copylist,
    const typename AT::tdual_float_2d &shearone,
    const int &dnum):
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _shearone(shearone.template view<DeviceType>()),
    _dnum(dnum)
  {
    _buf = typename AT::t_xfloat_1d_um(buf.template view<DeviceType>().data(),buf.extent(0)*buf.extent(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &mysend) const {
    const int i = _sendlist(mysend);
    int m = i*_dnum;
    for (int v = 0; v < _dnum; v++) {
      _buf(m++) = _shearone(i,v);
    }
    const int j = _copylist(mysend);
    if (j > -1) {
      for (int v = 0; v < _dnum; v++) {
	_shearone(i,v) = _shearone(j,v);
      }
    }
  }
 };

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int FixWallGranKokkos<DeviceType>::pack_exchange_kokkos(
  const int &nsend,
  DAT::tdual_xfloat_2d &buf,
  DAT::tdual_int_1d k_sendlist,
  DAT::tdual_int_1d k_copylist,
  ExecutionSpace space, int dim,
  X_FLOAT lo, X_FLOAT hi)
{
  k_shearone.template sync<DeviceType>();
  Kokkos::parallel_for(
    nsend,
    FixWallGranKokkos_PackExchangeFunctor<DeviceType>(
      buf,k_sendlist,k_copylist,k_shearone,sheardim));
  return nsend*sheardim;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct FixWallGranKokkos_UnpackExchangeFunctor
{
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_xfloat_1d_um _buf;
  typename AT::t_float_2d _shearone;
  typename AT::t_int_1d _indices;
  const int _dnum;

  FixWallGranKokkos_UnpackExchangeFunctor(
    const typename AT::tdual_xfloat_2d buf,
    const typename AT::tdual_float_2d &shearone,
    const typename AT::tdual_int_1d &indices,
    const int &dnum):
    _shearone(shearone.template view<DeviceType>()),
    _indices(indices.template view<DeviceType>()),
    _dnum(dnum)
  {
    _buf = typename AT::t_xfloat_1d_um(buf.template view<DeviceType>().data(),buf.extent(0)*buf.extent(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    int index = _indices(i);
    if (index > 0) {
      int m = i*_dnum;
      for (int v = 0; v < _dnum; v++) {
	_shearone(i,v) = _buf(m++);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf,
  DAT::tdual_int_1d &indices,int nrecv,
  int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
  ExecutionSpace space)
{
  Kokkos::parallel_for(
    nrecv/(atom->avec->size_border + atom->avec->size_velocity + 2),
    FixWallGranKokkos_UnpackExchangeFunctor<DeviceType>(
      k_buf,k_shearone,indices,sheardim));

  k_shearone.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixWallGranKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixWallGranKokkos<LMPHostType>;
#endif
}
