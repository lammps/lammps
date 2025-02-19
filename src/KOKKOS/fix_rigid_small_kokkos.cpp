// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_rigid_small_kokkos.h"

#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "atom_masks.h"
#include "molecule.h"
#include "math_const.h"
#include "math_eigen.h"
#include "math_extra.h"
#include "rigid_const.h"
#include "memory.h"
#include "update.h"
#include "force.h"
#include "random_mars.h"
#include "domain.h"
#include "error.h"
#include "math_extra_kokkos.h"
#include "kokkos_few.h"
#include "domain_kokkos.h"

#include <cmath>
#include <cstring>
#include <map>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

#define RVOUS 1   // 0 for irregular, 1 for all2all

void print_body(Body &b, int rank, int ibody, double **x){
  printf("==========\nrank %d, body %d\n", rank, ibody);
  printf("natoms = %d | ilocal = %d | xcm = %.6f %.6f %.6f\n", b.natoms, b.ilocal, b.xcm[0], b.xcm[1], b.xcm[2]);
  int i = b.ilocal;
  printf("owning position: %.6f %.6f %.6f\n", x[i][0], x[i][1], x[i][2]);
  printf("==========\n");
}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::FixRigidSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidSmall(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  commKK = (CommKokkos *) comm;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | V_MASK | VIRIAL_MASK | TYPE_MASK | TAG_MASK;
  datamask_modify = X_MASK | V_MASK | VIRIAL_MASK;

  maxexchange = 12 + bodysize;

  grow_arrays(atom->nmax);

  //TODO
  // for now, keep these 0 so communication uses super class methods
  //exchange_comm_device = 1;
  //forward_comm_device = 1;
  //reverse_comm_device = 1;
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::~FixRigidSmallKokkos()
{
  // TODO: Anything?
  if (copymode) return;

  atomKK->sync(Host, ALL_MASK);
  /*
  for(int ibody = 0; ibody < nlocal_body; ibody++){
    print_body(d_body(ibody), comm->me, ibody, atom->x);
  }
  fflush(stdout);
  */
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::init()
{
  FixRigidSmall::init();
  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot yet use respa with Kokkos");

  // TODO: Anything else?
}

/* ----------------------------------------------------------------------
   setup static/dynamic properties of rigid bodies, using current atom info.
   if reinitflag is not set, do the initialization only once, b/c properties
   may not be re-computable especially if overlapping particles or bodies
   are inserted from mol template.
     do not do dynamic init if read body properties from inpfile. this
   is b/c the inpfile defines the static and dynamic properties and may not
   be computable if contain overlapping particles setup_bodies_static()
   reads inpfile itself.
     cannot do this until now, b/c requires comm->setup() to have setup stencil
   invoke pre_neighbor() to ensure body xcmimage flags are reset
     needed if Verlet::setup::pbc() has remapped/migrated atoms for 2nd run
     setup_bodies_static() invokes pre_neighbor itself
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup_pre_neighbor()
{
  atomKK->sync(Host, datamask_read); // TODO: update to device

  FixRigidSmall::setup_pre_neighbor();

  atomKK->modified(Host, datamask_modify); // TODO: update to device
  atomKK->sync(execution_space, datamask_read); // TODO: update to device
}

/* ----------------------------------------------------------------------
   compute initial fcm and torque on bodies, also initial virial
   reset all particle velocities to be consistent with vcm and omega
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup(int vflag)
{
  atomKK->sync(Host, datamask_read); // TODO: update to device

  FixRigidSmall::setup(vflag);

  atomKK->modified(Host, datamask_modify); // TODO: update to device

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  auto h_bodyown = Kokkos::create_mirror_view(d_bodyown);
  auto h_bodytag = Kokkos::create_mirror_view(d_bodytag);
  auto h_atom2body = Kokkos::create_mirror_view(d_atom2body);
  auto h_xcmimage = Kokkos::create_mirror_view(d_xcmimage);
  auto h_displace = Kokkos::create_mirror_view(d_displace);
  auto h_vatom = Kokkos::create_mirror_view(d_vatom);
  for(int i = 0; i < nlocal+nghost; i++){ // TODO check local vs. ghost for all
    h_bodyown(i) = bodyown[i];
    h_bodytag(i) = bodytag[i];
    h_atom2body(i) = atom2body[i];
    h_xcmimage(i) = xcmimage[i];
    for(int j = 0; j < 3; j++){
      h_displace(i, j) = displace[i][j];
    }
    for(int j = 0; j < 6; j++){
      h_vatom(i, j) = vatom[i][j];
    }
  }
  Kokkos::deep_copy(d_bodyown, h_bodyown);
  Kokkos::deep_copy(d_bodytag, h_bodytag);
  Kokkos::deep_copy(d_atom2body, h_atom2body);
  Kokkos::deep_copy(d_xcmimage, h_xcmimage);
  Kokkos::deep_copy(d_displace, h_displace);
  Kokkos::deep_copy(d_vatom, h_vatom); //vatom unnecessary?

  Kokkos::resize(d_body, nmax_body);
  auto h_body = Kokkos::create_mirror_view(d_body);
  for(int i = 0; i < nlocal_body + nghost_body; i++){
    copy_body(&h_body(i), &body[i]);
  }
  Kokkos::deep_copy(d_body, h_body);

  // Before this, legacy communication routines were used
  // From now on, do device communication
  exchange_comm_device = 1;
  forward_comm_device = 1;
  sort_device = 1;
  reverse_comm_device = 1;
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::pre_neighbor(){
  if (!setupflag) {
    FixRigidSmall::pre_neighbor();
    return;
  }
  Kokkos::Profiling::pushRegion("rigid/small pre_neighbor");

  int triclinic = domain->triclinic;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  Few<double,6> h(domain->h);
  Few<double,6> h_inv(domain->h_inv);
  Few<double,3> boxlo(domain->boxlo);
  Few<double,3> lo, hi, period;

  if (triclinic == 0) {
    lo = Few<double,3>(domain->boxlo);
    hi = Few<double,3>(domain->boxhi);
    period = Few<double,3>(domain->prd);
  } else {
    lo = Few<double,3>(domain->boxlo_lamda);
    hi = Few<double,3>(domain->boxhi_lamda);
    period = Few<double,3>(domain->prd_lamda);
  }

  auto d_body = this->d_body;

  Kokkos::parallel_for(
    "fix rigid/small remap",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = d_body(ibody);
      Few<double,3> x(b.xcm);
      Few<double,3> coord;
      if (triclinic == 0) {
        coord = Few<double,3>(x);
      } else {
        coord = DomainKokkos::x2lamda(boxlo, h_inv, x);
      }
      // TODO: check
      DomainKokkos::remap(lo, hi, period, xperiodic, yperiodic, zperiodic, coord, b.image);
      //coord = DomainKokkos::remap(lo, hi, period, xperiodic, yperiodic, zperiodic, coord, b.image);
      if (triclinic) {
        x = DomainKokkos::lamda2x(boxlo, h, coord);
        for(int i = 0; i < 3; i++) {
          b.xcm[i] = x[i];
        }
      }
    }
  );

  nghost_body = 0;
  max_body_sent = 0;
  n_body_recv.clear();
  n_body_sent.clear();
  first_body.clear();
  commflag = BODY_SENDLIST;
  commKK->forward_comm_device<DeviceType>(this, 1);
  commflag = FULL_BODY;
  commKK->forward_comm_device<DeviceType>(this, bodysize);
  reset_atom2body();
  //check(4);

  int computed_nlocal_body = 0;
  int nlocal = atom->nlocal;
  auto d_bodyown = this->d_bodyown;
  int nlocal_body = this->nlocal_body;
  Kokkos::parallel_reduce(
    "pre_neighbor sanity check",
    Range1D(0, nlocal),
    KOKKOS_LAMBDA(const int i, int &count) {
      if (d_bodyown(i) >= 0) count++;
#ifndef KOKKOS_ENABLE_CUDA
      if (d_bodyown(i) >= 0 && nlocal_body==0) {
        error->one(FLERR, "atom {} has bodyown {} but no bodies", i, d_bodyown(i));
      }
      if(d_bodyown(i) >= nlocal_body){
        error->message(FLERR, "rank {} atom {} has bodyown {} but only {} local bodies",
            comm->me, i, d_bodyown(i), nlocal_body);
      }
#endif
    },
    computed_nlocal_body
  );
  if (nlocal_body != computed_nlocal_body) {
    printf("rank %d nlocal_body: %d vs %d\n", comm->me, nlocal_body, computed_nlocal_body);
    error->one(FLERR, "disagree!");
  }
  fflush(stdout);

  /*
  printf("rank %d bodytag:", comm->me);
  for(int i = 0; i < atom->nlocal; i++){
    if (d_bodytag(i) <= 0) continue;

    printf(" %d:%d", i, d_bodytag(i));
  }
  printf("\n");
  fflush(stdout);
  printf("rank %d atom2body:", comm->me);
  for(int i = 0; i < atom->nlocal; i++){
    if (d_atom2body(i) < 0) continue;

    printf(" %d:%d", i, d_atom2body(i));
  }
  printf("\n");
  fflush(stdout);
  */

  image_shift();
  Kokkos::Profiling::popRegion();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::initial_integrate(int vflag)
{
  Kokkos::Profiling::pushRegion("rigid/small initial integrate");
  copymode = 1;
  Kokkos::parallel_for("fix rigid/small/kk initial_integrate",
    Kokkos::RangePolicy<DeviceType, TagInitialIntegrate>(0, nlocal_body), *this
  );
  copymode = 0;

  // virial setup before call to set_xv

  v_init(vflag);

  // forward communicate updated info of all bodies

  commflag = INITIAL;
  commKK->forward_comm_device<DeviceType>(this,29);

  // set coords/orient and velocity/rotation of atoms in rigid bodies

  set_xv_kokkos(1);
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagInitialIntegrate, const int ibody) const {
  Body &b = d_body(ibody);

  // update vcm by 1/2 step

  double dtfm = dtf / b.mass;
  b.vcm[0] += dtfm * b.fcm[0];
  b.vcm[1] += dtfm * b.fcm[1];
  b.vcm[2] += dtfm * b.fcm[2];

  // update xcm by full step

  b.xcm[0] += dtv * b.vcm[0];
  b.xcm[1] += dtv * b.vcm[1];
  b.xcm[2] += dtv * b.vcm[2];

  // update angular momentum by 1/2 step

  b.angmom[0] += dtf * b.torque[0];
  b.angmom[1] += dtf * b.torque[1];
  b.angmom[2] += dtf * b.torque[2];

  // compute omega at 1/2 step from angmom at 1/2 step and current q
  // update quaternion a full step via Richardson iteration
  // returns new normalized quaternion, also updated omega at 1/2 step
  // update ex,ey,ez to reflect new quaternion

  MathExtraKokkos::angmom_to_omega(b.angmom,b.ex_space,b.ey_space,
                              b.ez_space,b.inertia,b.omega);
  MathExtraKokkos::richardson(b.quat,b.angmom,b.omega,b.inertia,dtq);
  MathExtraKokkos::q_to_exyz(b.quat,b.ex_space,b.ey_space,b.ez_space);
}

/* ----------------------------------------------------------------------
   apply Langevin thermostat to all 6 DOF of rigid bodies I own
   unlike fix langevin, this stores extra force in extra arrays,
     which are added in when a new fcm/torque are calculated
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::apply_langevin_thermostat_kokkos()
{
  copy_body_host();
  FixRigidSmall::apply_langevin_thermostat();
  auto h_langextra = Kokkos::create_mirror_view(d_langextra);
  for(int i = 0; i < nlocal_body; i++){
    for(int j = 0; j < 6; j++){
      h_langextra(i,j) = langextra[i][j];
    }
  }
  Kokkos::deep_copy(d_langextra, h_langextra);
}

/* ----------------------------------------------------------------------
   called from FixEnforce post_force() for 2d problems
   zero all body values that should be zero for 2d model
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::enforce2d()
{
  auto d_body = this->d_body;
  auto langflag = this->langflag && (langextra!=nullptr);
  auto d_langextra = this->d_langextra;

  Kokkos::parallel_for(
      "fix rigid/small zero z components",
      Range1D(0, nlocal_body),
      KOKKOS_LAMBDA (const int ibody) {
        Body &b = d_body(ibody);
        b.xcm[2] = 0.0;
        b.vcm[2] = 0.0;
        b.fcm[2] = 0.0;
        b.xgc[2] = 0.0;
        b.torque[0] = 0.0;
        b.torque[1] = 0.0;
        b.angmom[0] = 0.0;
        b.angmom[1] = 0.0;
        b.omega[0] = 0.0;
        b.omega[1] = 0.0;
        if(langflag) {
          d_langextra(ibody,2) = 0.0;
          d_langextra(ibody,3) = 0.0;
          d_langextra(ibody,4) = 0.0;
        }
      }
  );
  // TODO: Ensure unnecessary
  if (langflag) {
    for(int ibody = 0; ibody < nlocal_body; ibody++){
      langextra[ibody][2] = 0.0;
      langextra[ibody][3] = 0.0;
      langextra[ibody][4] = 0.0;
    }
  }
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::post_force(int /*vflag*/)
{
  Kokkos::Profiling::pushRegion("rigid/small post force");
  if (langflag) apply_langevin_thermostat_kokkos();
  if (earlyflag) compute_forces_and_torques_kokkos();
  Kokkos::Profiling::popRegion();
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::compute_forces_and_torques_kokkos()
{
  // sum over atoms to get force and torque on rigid body

  Kokkos::Profiling::pushRegion("rigix/small forces and torques");
  auto d_x = atomKK->k_x.view<DeviceType>();
  auto d_f = atomKK->k_f.view<DeviceType>();
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  auto d_body = this->d_body;

  Kokkos::parallel_for(
    "fix rigid/small zero fcm&tcm",
    Range1D(0,nlocal_body+nghost_body),
    KOKKOS_LAMBDA (const int ibody) {
      Body &b = d_body(ibody);
      b.fcm[0] = b.fcm[1] = b.fcm[2] = 0.0;
      b.torque[0] = b.torque[1] = b.torque[2] = 0.0;
    }
  );

  auto d_xcmimage = this->d_xcmimage;
  auto d_atom2body = this->d_atom2body;

  auto prd = Few<double,3>(domain->prd);
  auto h = Few<double,6>(domain->h);
  auto triclinic = domain->triclinic;

  Kokkos::parallel_for(
    "fix rigid/small add tcm&fcm",
    Range1D(0,nlocal),
    KOKKOS_LAMBDA (const int i) {
      if (d_atom2body(i) < 0) return;
      Body &b = d_body(d_atom2body(i));

      Kokkos::atomic_add(&b.fcm[0], d_f(i,0));
      Kokkos::atomic_add(&b.fcm[1], d_f(i,1));
      Kokkos::atomic_add(&b.fcm[2], d_f(i,2));

      Few<double,3> x_i;
      x_i[0] = d_x(i,0);
      x_i[1] = d_x(i,1);
      x_i[2] = d_x(i,2);
      auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,d_xcmimage(i));
      double dx = unwrap[0] - b.xcm[0];
      double dy = unwrap[1] - b.xcm[1];
      double dz = unwrap[2] - b.xcm[2];

      Kokkos::atomic_add(&b.torque[0], dy*d_f(i,2) - dz*d_f(i,1));
      Kokkos::atomic_add(&b.torque[1], dz*d_f(i,0) - dx*d_f(i,2));
      Kokkos::atomic_add(&b.torque[2], dx*d_f(i,1) - dy*d_f(i,0));
    }
  );


  // reverse communicate fcm, torque of all bodies

  Kokkos::Profiling::pushRegion("reverse communicate");
  commflag = FORCE_TORQUE;
  commKK->reverse_comm_device<DeviceType>(this,6);
  Kokkos::Profiling::popRegion();

  // include Langevin thermostat forces and torques

  // TODO: GPU langevin
  if (langflag) {
    Kokkos::Profiling::pushRegion("rigid/small langevin");
    copy_body_host();
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      double *fcm = body[ibody].fcm;
      fcm[0] += langextra[ibody][0];
      fcm[1] += langextra[ibody][1];
      fcm[2] += langextra[ibody][2];
      double *tcm = body[ibody].torque;
      tcm[0] += langextra[ibody][3];
      tcm[1] += langextra[ibody][4];
      tcm[2] += langextra[ibody][5];
    }
    copy_body_device();
    Kokkos::Profiling::popRegion();
  }

  // add gravity force to COM of each body

  // TODO
  if (id_gravity) {
    error->all(FLERR, "gravity not implemented");
    double mass;
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      mass = body[ibody].mass;
      double *fcm = body[ibody].fcm;
      fcm[0] += gvec[0]*mass;
      fcm[1] += gvec[1]*mass;
      fcm[2] += gvec[2]*mass;
    }
  }

  Kokkos::Profiling::popRegion();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::final_integrate()
{
  Kokkos::Profiling::pushRegion("rigid/small final integrate");

  //check(3);

  if (!earlyflag) compute_forces_and_torques_kokkos();

  // update vcm and angmom, recompute omega

  auto nlocal_body = this->nlocal_body;
  auto d_body = this->d_body;
  auto dtf = this->dtf;

  Kokkos::parallel_for(
    "fix rigid/small update vcm+angmom+omega",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA (const int ibody) {
      Body &b = d_body(ibody);

      // update vcm by 1/2 step

      double dtfm = dtf / b.mass;
      b.vcm[0] += dtfm * b.fcm[0];
      b.vcm[1] += dtfm * b.fcm[1];
      b.vcm[2] += dtfm * b.fcm[2];

      // update angular momentum by 1/2 step

      b.angmom[0] += dtf * b.torque[0];
      b.angmom[1] += dtf * b.torque[1];
      b.angmom[2] += dtf * b.torque[2];

      MathExtraKokkos::angmom_to_omega(&b.angmom[0],&b.ex_space[0],&b.ey_space[0],
                                &b.ez_space[0],&b.inertia[0],&b.omega[0]);
    }
  );

  // forward communicate updated info of all bodies

  commflag = FINAL;
  commKK->forward_comm_device<DeviceType>(this,10);

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_xv_kokkos(0);
  Kokkos::Profiling::popRegion();
}

/* ----------------------------------------------------------------------
   count # of DOF removed by rigid bodies for atoms in igroup
   return total count of DOF
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::dof(int tgroup)
{
  if (!setupflag) {
    int nlocal = atom->nlocal;
    copy_body_host();
    auto h_bodyown = Kokkos::create_mirror_view_and_copy(LMPHostType(), d_bodyown);
    auto h_bodytag = Kokkos::create_mirror_view_and_copy(LMPHostType(), d_bodytag);
    memcpy(bodyown, h_bodyown.data(), nlocal*sizeof(int));
    memcpy(bodytag, h_bodytag.data(), nlocal*sizeof(int));
  }

  return FixRigidSmall::dof(tgroup);
}

/* ----------------------------------------------------------------------
   adjust xcm of each rigid body due to box deformation
   called by various fixes that change box size/shape
   flag = 0/1 means map from box to lamda coords or vice versa
------------------------------------------------------------------------- */

// TODO: deviceify
template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::deform(int flag)
{
  copy_body_host();
  if (flag == 0)
    for (int ibody = 0; ibody < nlocal_body; ibody++)
      domain->x2lamda(body[ibody].xcm,body[ibody].xcm);
  else
    for (int ibody = 0; ibody < nlocal_body; ibody++)
      domain->lamda2x(body[ibody].xcm,body[ibody].xcm);
  copy_body_device();
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
   setxflag = whether to update positions
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_xv_kokkos(int setxflag)
{
  Kokkos::Profiling::pushRegion("rigid/small set_xv");
  this->xprd = domain->xprd;
  this->yprd = domain->yprd;
  this->zprd = domain->zprd;
  this->xy = domain->xy;
  this->xz = domain->xz;
  this->yz = domain->yz;

  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  int nlocal = atom->nlocal;

  EV_FLOAT ev;
  if(vflag_atom){
    Kokkos::deep_copy(d_vatom, 0.0);
  }
  copymode = 1;
  if (setxflag) {
    Kokkos::parallel_reduce("fix rigid/small setxv", Kokkos::RangePolicy<DeviceType, TagSetXV<1>>(0, nlocal), *this, ev);
    Kokkos::parallel_for("fix rigid/small update xgc", Kokkos::RangePolicy<DeviceType, TagUpdateXGC>(0, nlocal_body+nghost_body), *this);
  } else {
    Kokkos::parallel_reduce("fix rigid/small setv", Kokkos::RangePolicy<DeviceType, TagSetXV<0>>(0, nlocal), *this, ev);
  }
  copymode = 0;
  if(evflag){
    if(vflag_global){
      for(int i = 0; i < 6; i++){
        virial[i] += ev.v[i];
      }
    }
    if(vflag_atom){
      auto h_vatom = Kokkos::create_mirror_view_and_copy(LMPHostType(), d_vatom);
      for(int i = 0; i < nlocal; i++){
        for(int j = 0; j < 6; j++){
          vatom[i][j] = h_vatom(i, j);
        }
      }
    }
  }
  // TODO: Specialize
  atomKK->modified(execution_space, datamask_modify);
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
template<int SETXFLAG>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagSetXV<SETXFLAG>, const int i, EV_FLOAT &ev) const{
  if (d_atom2body(i) < 0) return;
  Body &b = d_body(d_atom2body(i));

  double xbox = (d_xcmimage(i) & IMGMASK) - IMGMAX;
  double ybox = (d_xcmimage(i) >> IMGBITS & IMGMASK) - IMGMAX;
  double zbox = (d_xcmimage(i) >> IMG2BITS) - IMGMAX;

  double x0, x1, x2, v0, v1, v2;

  // save old positions and velocities for virial

  if (evflag) {
    if (triclinic == 0) {
      x0 = d_x(i,0) + xbox*xprd;
      x1 = d_x(i,1) + ybox*yprd;
      x2 = d_x(i,2) + zbox*zprd;
    } else {
      x0 = d_x(i,0) + xbox*xprd + ybox*xy + zbox*xz;
      x1 = d_x(i,1) + ybox*yprd + zbox*yz;
      x2 = d_x(i,2) + zbox*zprd;
    }
    v0 = d_v(i,0);
    v1 = d_v(i,1);
    v2 = d_v(i,2);
  }

  // x = displacement from center-of-mass, based on body orientation
  // v = vcm + omega around center-of-mass

  double delta[3];
  MathExtraKokkos::matvec(b.ex_space,b.ey_space,b.ez_space,&d_displace(i,0),delta);

  d_v(i,0) = b.omega[1]*delta[2] - b.omega[2]*delta[1] + b.vcm[0];
  d_v(i,1) = b.omega[2]*delta[0] - b.omega[0]*delta[2] + b.vcm[1];
  d_v(i,2) = b.omega[0]*delta[1] - b.omega[1]*delta[0] + b.vcm[2];

  // add center of mass to displacement
  // map back into periodic box via xbox,ybox,zbox
  // for triclinic, add in box tilt factors as well

  if constexpr(SETXFLAG) {
    if (triclinic == 0) {
      d_x(i,0) = delta[0] + b.xcm[0] - xbox*xprd;
      d_x(i,1) = delta[1] + b.xcm[1] - ybox*yprd;
      d_x(i,2) = delta[2] + b.xcm[2] - zbox*zprd;
    } else {
      d_x(i,0) = delta[0] + b.xcm[0] - xbox*xprd - ybox*xy - zbox*xz;
      d_x(i,1) = delta[1] + b.xcm[1] - ybox*yprd - zbox*yz;
      d_x(i,2) = delta[2] + b.xcm[2] - zbox*zprd;
    }
  }

  // virial = unwrapped coords dotted into body constraint force
  // body constraint force = implied force due to v change minus f external
  // assume f does not include forces internal to body
  // 1/2 factor b/c final_integrate contributes other half
  // assume per-atom contribution is due to constraint force on that atom

  if (evflag) { // TODO: Figure out
    double massone;
    if (d_rmass.data()) massone = d_rmass(i);
    else massone = d_mass(d_type(i));
    double fc0 = massone*(d_v(i,0) - v0)/dtf - d_f(i,0);
    double fc1 = massone*(d_v(i,1) - v1)/dtf - d_f(i,1);
    double fc2 = massone*(d_v(i,2) - v2)/dtf - d_f(i,2);

    double vr[6];
    vr[0] = 0.5*x0*fc0;
    vr[1] = 0.5*x1*fc1;
    vr[2] = 0.5*x2*fc2;
    vr[3] = 0.5*x0*fc1;
    vr[4] = 0.5*x0*fc2;
    vr[5] = 0.5*x1*fc2;

    double rlist[3] = {x0, x1, x2};
    double flist[3] = {0.5*fc0, 0.5*fc1, 0.5*fc2};
    v_tally(ev,i,vr,rlist,flist,b.xgc);
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::operator()(TagUpdateXGC, const int ibody) const {
  Body &b = d_body(ibody);
  MathExtraKokkos::matvec(b.ex_space,b.ey_space,b.ez_space,
                    b.xgc_body,b.xgc);
  b.xgc[0] += b.xcm[0];
  b.xgc[1] += b.xcm[1];
  b.xgc[2] += b.xcm[2];
}


/* ----------------------------------------------------------------------
   write out restart info for mass, COM, inertia tensor to file
   identical format to inpfile option, so info can be read in when restarting
   each proc contributes info for rigid bodies it owns
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::write_restart_file(const char *file)
{
  copy_body_host();
  FixRigidSmall::write_restart_file(file);
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_arrays(int nmax)
{
  Kokkos::Profiling::pushRegion("rigid/small grow arrays");
  int prev_size = d_bodyown.extent(0);
  FixRigidSmall::grow_arrays(nmax);
  Kokkos::resize(d_bodyown, nmax);
  Kokkos::resize(d_bodytag, nmax);
  Kokkos::resize(d_atom2body, nmax);
  Kokkos::resize(d_xcmimage, nmax);
  Kokkos::resize(d_displace, nmax, 3);
  Kokkos::resize(d_vatom, nmax, 6);

  auto d_bodyown = this->d_bodyown;
  auto d_atom2body = this->d_atom2body;
  Kokkos::parallel_for(
    "fix rigid/small set new bodyown and atom2body to -1",
    Range1D(prev_size, nmax),
    KOKKOS_LAMBDA(const int i) {
      d_bodyown(i) = -1;
      d_atom2body(i) = -1;
    }
  );
  Kokkos::Profiling::popRegion();
}

/* ----------------------------------------------------------------------
   initialize a molecule inserted by another fix, e.g. deposit or pour
   called when molecule is created
   nlocalprev = # of atoms on this proc before molecule inserted
   tagprev = atom ID previous to new atoms in the molecule
   xgeom = geometric center of new molecule
   vcm = COM velocity of new molecule
   quat = rotation of new molecule (around geometric center)
          relative to template in Molecule class
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_molecule(int nlocalprev, tagint tagprev, int imol,
                                 double *xgeom, double *vcm, double *quat)
{
  // TODO: Eliminate host work

  // copy current bodies to host because they will later overwrite device
  copy_body_host();

  // add new bodies on host
  FixRigidSmall::set_molecule(nlocalprev, tagprev, imol, xgeom, vcm, quat);
  // update device body list to same size and copy
  grow_body();
  copy_body_device();

  int nlocal = atom->nlocal;

  // copy bodytag/own of new atoms
  auto h_bodytag = Kokkos::create_mirror_view(d_bodytag);
  auto h_bodyown = Kokkos::create_mirror_view(d_bodyown);

  for(int i = nlocalprev; i < nlocal; i++){
    h_bodytag(i) = bodytag[i];
    h_bodyown(i) = bodyown[i];
  }
  Kokkos::deep_copy(d_bodytag, h_bodytag);
  Kokkos::deep_copy(d_bodyown, h_bodyown);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_exchange_kokkos (
    const int &nsend,DAT::tdual_xfloat_2d &k_buf,
    DAT::tdual_int_1d k_sendlist,
    DAT::tdual_int_1d k_copylist,
    ExecutionSpace space) {

  Kokkos::Profiling::pushRegion("rigid/small pack exchange");

  k_buf.sync<DeviceType>();
  k_copylist.sync<DeviceType>();
  k_sendlist.sync<DeviceType>();

  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_copylist = k_copylist.view<DeviceType>();
  auto d_sendlist = k_sendlist.view<DeviceType>();

  auto d_vatom = this->d_vatom;
  auto d_bodytag = this->d_bodytag;
  auto d_xcmimage = this->d_xcmimage;
  auto d_displace = this->d_displace;
  auto d_atom2body = this->d_atom2body;
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;
  auto exchange_size = this->maxexchange;

  int n_deleted_bodies = 0;

  // TODO: Optimize with parallel scan,
  // see fix shake/kk
  Kokkos::parallel_reduce(
    "fix rigid/small pack exchange",
    Range1D(0, nsend),
    KOKKOS_LAMBDA(const int isend, int &n_deleted_bodies) {
      const int i = d_sendlist(isend);
      int m = isend*exchange_size;

      d_buf(m++) = d_ubuf(d_bodytag(i)).d;
      if (d_bodytag(i)){
        d_buf(m++) = d_ubuf(d_xcmimage(i)).d;
        d_buf(m++) = d_displace(i,0);
        d_buf(m++) = d_displace(i,1);
        d_buf(m++) = d_displace(i,2);
        for (int k = 0; k < 6; k++) {
          d_buf(m++) = d_vatom(i,k);
        }

        if (d_bodyown(i) < 0) {
          d_buf(m++) = 0;
        }
        else {
          d_buf(m++) = 1;
          memcpy(&d_buf(m),&d_body(d_bodyown(i)),sizeof(Body));
        }
        if(d_bodyown(i) >= 0){
          // mark this slot as free for incoming bodies
          d_body(d_bodyown(i)).natoms = -1;
          n_deleted_bodies++;
          // d_bodyown(i) = -1;
          // d_bodytag(i) = -1;
          // d_atom2body(i) = -1;
        }
      }

      const int j = d_copylist(isend);
      if (j < 0) return;

      d_bodytag(i) = d_bodytag(j);
      d_xcmimage(i) = d_xcmimage(j);
      d_bodyown(i) = d_bodyown(j);
#ifndef KOKKOS_ENABLE_CUDA
      if(d_bodyown(i) >= nlocal_body){
        error->one(FLERR, "rank {} atom {} has bodyown {} but nlocal body {}",
            comm->me, i, d_bodyown(i), nlocal_body);
      }
#endif
      for(int k = 0; k < 3; k++)
        d_displace(i,k) = d_displace(j,k);
      for(int k = 0; k < 6; k++)
        d_vatom(i,k) = d_vatom(j,k);

      if (d_bodyown(i) >= 0) {
        d_body(d_bodyown(i)).ilocal = i;
      }

      // this appears necessary
      d_bodyown(j) = -1;
      d_bodytag(j) = -1;
      d_atom2body(j) = -1;
    },
    n_deleted_bodies
  );

  // Need to pack remaining bodies densely
  int new_nlocal_body = nlocal_body - n_deleted_bodies;
  IntView1D from_indices("from idx", n_deleted_bodies);

  // count bodies that need to be moved
  // and do cumulative sum to determine
  // their relative position
  int n_to_move = 0;
  Kokkos::parallel_scan(
    "fix rigid/small count bodies to move",
    Range1D(new_nlocal_body, nlocal_body),
    KOKKOS_LAMBDA(const int ibody, int &count, const bool is_final){
      if (is_final && d_body(ibody).natoms > 0) {
        from_indices(count) = ibody;
        // error->message(FLERR, "rank {} needs to move body from {}, this is #{}", comm->me, ibody, count);
      }
      if (d_body(ibody).natoms > 0) count++;
    },
    n_to_move
  );

  // printf("rank %d has %d bodies, deleting %d, moving %d\n", comm->me, nlocal_body, n_deleted_bodies, n_to_move);
  fflush(stdout);


  // count open slots, fill in corresponding body
  // need copylist for updating bodyown
  Kokkos::parallel_scan(
    "fix rigid/small create body copylist",
    Range1D(0, new_nlocal_body),
    KOKKOS_LAMBDA(const int i, int &count, const bool is_final){
      if (d_body(i).natoms<0) {
        if (is_final) {
          copy_body(&d_body(i), &d_body(from_indices(count)));
          d_bodyown(d_body(i).ilocal) = i;
          // error->message(FLERR, "rank {} moving body #{} from {} to {}, nlocal={}, ndeleted={}, ntomove={}", comm->me, count, from_indices(count), i, nlocal_body, n_deleted_bodies, n_to_move);
          // fflush(stdout);
        }
        count++;
      }
    }
  );

  nlocal_body = new_nlocal_body;
  Kokkos::Profiling::popRegion();

  return nsend*exchange_size;
}


/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &k_indices,int nrecv,
                              ExecutionSpace space)
{
  Kokkos::Profiling::pushRegion("rigid/small unpack exchange");
  k_buf.sync<DeviceType>();
  k_indices.sync<DeviceType>();

  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_indices = k_indices.view<DeviceType>();

  auto d_bodytag = this->d_bodytag;
  auto d_xcmimage = this->d_xcmimage;
  auto d_displace = this->d_displace;
  auto d_vatom = this->d_vatom;
  auto d_bodyown = this->d_bodyown;
  int exchange_size = this->maxexchange;

  // TODO: make this parallel scan to get body indices,
  // then do parallel_for over new bodies
  int n_new_body = 0;
  Kokkos::parallel_reduce(
    "fix rigid/small count incoming bodies",
    Range1D(0, nrecv),
    KOKKOS_LAMBDA(const int irecv, int &count){
      int m = irecv*exchange_size;
      int i = d_indices(irecv);
      if (i < 0) return;

      d_bodyown(i) = -1; // set later
      d_bodytag(i) = (tagint) d_ubuf(d_buf(m++)).i;
      if (d_bodytag(i)) {
        d_xcmimage(i) = (imageint) d_ubuf(d_buf(m++)).i;
        for(int k = 0; k < 3; k++){
          d_displace(i,k) = d_buf(m++);
        }
        for(int k = 0; k < 6; k++){
          d_vatom(i,k) = d_buf(m++);
        }
        if (d_buf(m++) > 0) count++;
      }
    },
    n_new_body
  );

  // printf("rank %d has %d bodies, receiving %d\n", comm->me, nlocal_body, n_new_body);
  // fflush(stdout);

  while (nlocal_body + n_new_body > nmax_body) {
    grow_body();
  }
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;

  Kokkos::parallel_scan(
    "fix rigid/small copy incoming bodies",
    Range1D(0, nrecv),
    KOKKOS_LAMBDA(const int irecv, int &count, const bool is_final){
      int i = d_indices(irecv);
      if(i < 0) return;

      if(d_bodytag(i) <= 0) return;

      int m = irecv*exchange_size + 11;

      // if owning body
      if (d_buf(m) > 0) {
        if (is_final) {
          int body_idx = nlocal_body + count;
          d_bodyown(i) = body_idx;
          memcpy(&d_body(body_idx),&d_buf(m+1),sizeof(Body));
          d_body(body_idx).ilocal = i;
        }
        count++;
      }
    }
  );


  // error->warning(FLERR, "rank {} receiving {} bodies", comm->me, n_new_body);
  this->nlocal_body += n_new_body;

  /*
  printf("rank %d now has %d local bodies\n", comm->me, nlocal_body);
  for(int i = 0; i < nlocal_body; i++){
    print_body(d_body(i), comm->me, i, atom->x);
  }
  fflush(stdout);
  */
  Kokkos::Profiling::popRegion();
}

/* ----------------------------------------------------------------------
   only pack body info if own or ghost atom owns the body
   for FULL_BODY, send 0/1 flag with every atom
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                                                        int iswap, DAT::tdual_xfloat_1d &k_buf,
                                                        int pbc_flag, int* pbc)

{
  Kokkos::Profiling::pushRegion("rigid/small pack forward");

  auto d_sendlist = Kokkos::subview(k_sendlist.view<DeviceType>(), iswap, Kokkos::ALL);
  auto d_buf = k_buf.view<DeviceType>();

  int m = 0;

  auto d_body = this->d_body;
  auto d_bodyown = this->d_bodyown;

  if (commflag == INITIAL) {
    auto bodysize = this->bodysize;
    int n_body = n_body_sent[iswap];
    auto d_body_sendlist = d_body_sendlists[iswap];
    Kokkos::parallel_for("fix rigid/small pack forward comm initial",
      Range1D(0, n_body),
      KOKKOS_LAMBDA (const int ibodysend) {
        int ibody = d_body_sendlist(ibodysend);
        Body &b = d_body(ibody);
        int m = 29*ibodysend;
        d_buf(m++) = b.xcm[0];
        d_buf(m++) = b.xcm[1];
        d_buf(m++) = b.xcm[2];
        d_buf(m++) = b.xgc[0];
        d_buf(m++) = b.xgc[1];
        d_buf(m++) = b.xgc[2];
        d_buf(m++) = b.vcm[0];
        d_buf(m++) = b.vcm[1];
        d_buf(m++) = b.vcm[2];
        d_buf(m++) = b.quat[0];
        d_buf(m++) = b.quat[1];
        d_buf(m++) = b.quat[2];
        d_buf(m++) = b.quat[3];
        d_buf(m++) = b.omega[0];
        d_buf(m++) = b.omega[1];
        d_buf(m++) = b.omega[2];
        d_buf(m++) = b.ex_space[0];
        d_buf(m++) = b.ex_space[1];
        d_buf(m++) = b.ex_space[2];
        d_buf(m++) = b.ey_space[0];
        d_buf(m++) = b.ey_space[1];
        d_buf(m++) = b.ey_space[2];
        d_buf(m++) = b.ez_space[0];
        d_buf(m++) = b.ez_space[1];
        d_buf(m++) = b.ez_space[2];
        d_buf(m++) = b.conjqm[0];
        d_buf(m++) = b.conjqm[1];
        d_buf(m++) = b.conjqm[2];
        d_buf(m++) = b.conjqm[3];
      }
    );
    Kokkos::Profiling::popRegion();
    return 29*n_body;
  }
  if (commflag == FINAL) {
    auto bodysize = this->bodysize;
    int n_body = n_body_sent[iswap];
    auto d_body_sendlist = d_body_sendlists[iswap];
    Kokkos::parallel_for("fix rigid/small pack forward comm final",
      Range1D(0, n_body),
      KOKKOS_LAMBDA (const int ibodysend) {
        int ibody = d_body_sendlist(ibodysend);
        Body &b = d_body(ibody);
        int m = 10*ibodysend;
        d_buf(m++) = b.vcm[0];
        d_buf(m++) = b.vcm[1];
        d_buf(m++) = b.vcm[2];
        d_buf(m++) = b.omega[0];
        d_buf(m++) = b.omega[1];
        d_buf(m++) = b.omega[2];
        d_buf(m++) = b.conjqm[0];
        d_buf(m++) = b.conjqm[1];
        d_buf(m++) = b.conjqm[2];
        d_buf(m++) = b.conjqm[3];
      }
    );
    Kokkos::Profiling::popRegion();
    return 10*n_body;
  } else if (commflag == FULL_BODY) {
    auto bodysize = this->bodysize;
    int n_body = n_body_sent[iswap];
    auto d_body_sendlist = d_body_sendlists[iswap];
    Kokkos::parallel_for(
      "fix rigid/small full body pack forward",
      Range1D(0, n_body),
      KOKKOS_LAMBDA(const int ibodysend) {
        int ibody = d_body_sendlist(ibodysend);
        int m = (bodysize)*ibodysend;

        memcpy(&d_buf(m), &d_body(ibody), sizeof(Body));
      }
    );
    Kokkos::Profiling::popRegion();
    return bodysize*n_body;
  } else if (commflag == BODY_SENDLIST) {
    int n_sent = 0;
    // TODO: parallel_scan
    Kokkos::parallel_reduce(
      "fix rigid/small count bodies sent",
      Range1D(0, n),
      KOKKOS_LAMBDA(const int isend, int &tmp) {
        int i = d_sendlist(isend);
        if (d_bodyown(i) >= 0) {
          tmp++;
          d_buf(isend) = 1;
        } else {
          d_buf(isend) = 0;
        }
      },
      n_sent
    );
    n_body_sent[iswap] = n_sent;
    if (n_sent > max_body_sent) max_body_sent = n_sent;

    if (d_body_sendlists.count(iswap)==0 || d_body_sendlists[iswap].extent_int(0)<n_sent) {
      d_body_sendlists[iswap] = IntView1D("body sendlist", n_sent);
    }
    auto d_body_sendlist = d_body_sendlists[iswap];

    Kokkos::parallel_scan(
      "fix rigid/small create body sendlist",
      Range1D(0, n),
      KOKKOS_LAMBDA(const int isend, int &count, const bool is_final) {
        const int i = d_sendlist(isend);
        if (d_bodyown(i) >= 0) {
          if (is_final) {
            d_body_sendlist(count) = d_bodyown(i);
          }
          count++;
        }
      }
    );
    Kokkos::Profiling::popRegion();
    return n;
  }

  Kokkos::Profiling::popRegion();
  return m;
}

/* ----------------------------------------------------------------------
   only ghost atoms are looped over
   for FULL_BODY, store a new ghost body if this atom owns it
   for other commflag values, only unpack body info if atom owns it
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first, DAT::tdual_xfloat_1d &k_buf)
{
  if (n==0) return;
  Kokkos::Profiling::pushRegion("rigid/small unpack forward");
  this->first = first;
  auto d_buf = k_buf.view<DeviceType>();
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;

  if (commflag == INITIAL) {
    int n_incoming_bodies = n_body_recv[first];
    int start_body = first_body[first];
    Kokkos::parallel_for("fix rigid/small unpack forward comm initial",
      Range1D(0, n_incoming_bodies),
      KOKKOS_LAMBDA(const int ibodyrecv) {
        int ibody = ibodyrecv + start_body;
        int m = 29*ibodyrecv;
        Body &b = d_body(ibody);
        b.xcm[0] = d_buf(m++);
        b.xcm[1] = d_buf(m++);
        b.xcm[2] = d_buf(m++);
        b.xgc[0] = d_buf(m++);
        b.xgc[1] = d_buf(m++);
        b.xgc[2] = d_buf(m++);
        b.vcm[0] = d_buf(m++);
        b.vcm[1] = d_buf(m++);
        b.vcm[2] = d_buf(m++);
        b.quat[0] = d_buf(m++);
        b.quat[1] = d_buf(m++);
        b.quat[2] = d_buf(m++);
        b.quat[3] = d_buf(m++);
        b.omega[0] = d_buf(m++);
        b.omega[1] = d_buf(m++);
        b.omega[2] = d_buf(m++);
        b.ex_space[0] = d_buf(m++);
        b.ex_space[1] = d_buf(m++);
        b.ex_space[2] = d_buf(m++);
        b.ey_space[0] = d_buf(m++);
        b.ey_space[1] = d_buf(m++);
        b.ey_space[2] = d_buf(m++);
        b.ez_space[0] = d_buf(m++);
        b.ez_space[1] = d_buf(m++);
        b.ez_space[2] = d_buf(m++);
        b.conjqm[0] = d_buf(m++);
        b.conjqm[1] = d_buf(m++);
        b.conjqm[2] = d_buf(m++);
        b.conjqm[3] = d_buf(m++);
      }
    );
  } else if (commflag == FINAL) {
    int n_incoming_bodies = n_body_recv[first];
    int start_body = first_body[first];
    Kokkos::parallel_for("fix rigid/small/kk unpack forward comm final",
      Range1D(0, n_incoming_bodies),
      KOKKOS_LAMBDA(const int ibodyrecv) {
        int ibody = ibodyrecv + start_body;
        int m = 10*ibodyrecv;
        Body &b = d_body(ibody);
        b.vcm[0] = d_buf(m++);
        b.vcm[1] = d_buf(m++);
        b.vcm[2] = d_buf(m++);
        b.omega[0] = d_buf(m++);
        b.omega[1] = d_buf(m++);
        b.omega[2] = d_buf(m++);
        b.conjqm[0] = d_buf(m++);
        b.conjqm[1] = d_buf(m++);
        b.conjqm[2] = d_buf(m++);
        b.conjqm[3] = d_buf(m++);
      }
    );
  } else if (commflag == FULL_BODY) {
    Kokkos::Profiling::pushRegion("unpack forward full body");
    auto bodysize = this->bodysize;
    int n_incoming_bodies = n_body_recv[first];
    int start_body = first_body[first];

    Kokkos::parallel_for(
      "fix rigid/small pack incoming ghost bodies",
      Range1D(0, n_incoming_bodies),
      KOKKOS_LAMBDA(const int ibodyrecv) {
        int m = ibodyrecv*(bodysize);
        memcpy(&d_body(ibodyrecv+start_body), &d_buf(m), sizeof(Body));
      }
    );
    Kokkos::Profiling::popRegion();
  } else if (commflag == BODY_SENDLIST) {
    Kokkos::Profiling::pushRegion("unpack forward body sendlist");
    auto bodysize = this->bodysize;
    int n_curr_bodies = this->nlocal_body + this->nghost_body;
    first_body[first] = n_curr_bodies;
    int n_incoming_bodies = 0;
    Kokkos::parallel_scan(
      "fix rigid/small count incoming bodies",
      Range1D(0, n),
      KOKKOS_LAMBDA(const int irecv, int &count, const bool is_final) {
        int i = irecv+first;
        int m = irecv;
        d_bodyown(i) = d_buf(m);
        if (d_bodyown(i)) {
          if (is_final) d_bodyown(i) = n_curr_bodies + count;
          count++;
        } else {
          d_bodyown(i) = -1;
        }
      },
      n_incoming_bodies
    );
    if (n_body_recv.count(first)) {
      error->one(FLERR, "first={} should not already be key, receiving {} atoms {} bodies", first, n, n_incoming_bodies);
    }
    n_body_recv[first] = n_incoming_bodies;
    while (n_curr_bodies+n_incoming_bodies > nmax_body) {
      grow_body();
    }
    d_body = this->d_body;

    this->nghost_body += n_incoming_bodies;
    Kokkos::Profiling::popRegion();
  }
  Kokkos::Profiling::popRegion();
}
template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first, DAT::tdual_xfloat_1d &k_buf)
{
  if (commflag != FORCE_TORQUE) {
    error->all(FLERR, "attempting invalid reverse comm on device");
  }

  if (n==0) return 0;

  auto d_buf = k_buf.view<DeviceType>();
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;

  int n_body = n_body_recv[first];
  int start_body = first_body[first];

  Kokkos::parallel_for(
    "fix rigid/small pack reverse comm",
    Range1D(0, n_body),
    KOKKOS_LAMBDA(const int ibodysend) {
      int ibody = ibodysend + start_body;
      int m = ibodysend*6;
      Body &b = d_body(ibody);

      d_buf(m++) = b.fcm[0];
      d_buf(m++) = b.fcm[1];
      d_buf(m++) = b.fcm[2];
      d_buf(m++) = b.torque[0];
      d_buf(m++) = b.torque[1];
      d_buf(m++) = b.torque[2];
    }
  );

  return 6*n_body;
}


template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                                          int iswap, DAT::tdual_xfloat_1d &k_buf)
{
  if (commflag != FORCE_TORQUE) {
    error->all(FLERR, "attempting invalid reverse comm on device");
  }

  auto d_buf = k_buf.view<DeviceType>();
  auto d_bodyown = this->d_bodyown;
  auto d_body = this->d_body;
  auto d_sendlist = Kokkos::subview(k_sendlist.view<DeviceType>(), iswap, Kokkos::ALL);

  int n_body = n_body_sent[iswap];
  auto d_body_sendlist = d_body_sendlists[iswap];
  Kokkos::parallel_for(
    "fix rigid/small unpack reverse comm",
    Range1D(0, n_body),
    KOKKOS_LAMBDA (const int ibodyrecv) {
      int ibody = d_body_sendlist(ibodyrecv);
      Body &b = d_body(ibody);
      int m = 6*ibodyrecv;

      b.fcm[0] += d_buf(m++);
      b.fcm[1] += d_buf(m++);
      b.fcm[2] += d_buf(m++);
      b.torque[0] += d_buf(m++);
      b.torque[1] += d_buf(m++);
      b.torque[2] += d_buf(m++);
    }
  );
}

/* ----------------------------------------------------------------------
   grow body data structure
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_body()
{
  Kokkos::Profiling::pushRegion("rigid/small grow body");

  // In set_molecule, CPU code calls grow_body first
  // Want to not double grow
  // Still need to do it during initial FULL_BODY unpack forward comm
  if (nmax_body == d_body.extent_int(0) || d_body.extent_int(0) == 0) {
    FixRigidSmall::grow_body();
  }
  Kokkos::resize(d_body, nmax_body);
  if (langflag) Kokkos::resize(d_langextra, nmax_body,6);
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::reset_atom2body()
{
  if (!setupflag) {
    // called during setup_bodies
    FixRigidSmall::reset_atom2body();
    return;
  }
  Kokkos::Profiling::pushRegion("rigid/small reset_atom2body");

  int nlocal = atom->nlocal;

  auto d_atom2body = this->d_atom2body;
  auto d_bodytag = this->d_bodytag;
  auto d_bodyown = this->d_bodyown;
  auto atomKK = this->atomKK;

  auto map_style = atom->map_style;
  decltype(atomKK->k_map_array) k_map_array;
  decltype(atomKK->k_map_hash) k_map_hash;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }


  Kokkos::parallel_for(
    "fix rigid/small reset atom2body",
    Range1D(0, nlocal),
    KOKKOS_LAMBDA(const int i){
      d_atom2body(i) = -1;
      if (d_bodytag(i)) {
        // int iowner = atomKK->map(d_bodytag(i));
        int iowner = AtomKokkos::map_kokkos<DeviceType>(d_bodytag(i),map_style,k_map_array,k_map_hash);
        d_atom2body(i) = d_bodyown(iowner);
      }
    }
  );
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::image_shift()
{
  if (!setupflag) {
    // called during setup_bodies
    FixRigidSmall::image_shift();
    return;
  }
  Kokkos::Profiling::pushRegion("rigid/small image shift");

  ImageIntView1D d_image = atomKK->k_image.view<DeviceType>();
  int nlocal = atom->nlocal;
  auto d_body = this->d_body;
  auto d_xcmimage = this->d_xcmimage;
  auto d_atom2body = this->d_atom2body;

  Kokkos::parallel_for(
    "fix rigid/small image shift",
    Range1D(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      if (d_atom2body(i) < 0) return;
      Body &b = d_body(d_atom2body(i));

      imageint tdim,bdim,xdim[3];
      tdim = d_image(i) & IMGMASK;
      bdim = b.image & IMGMASK;
      xdim[0] = IMGMAX + tdim - bdim;
      tdim = (d_image(i) >> IMGBITS) & IMGMASK;
      bdim = (b.image >> IMGBITS) & IMGMASK;
      xdim[1] = IMGMAX + tdim - bdim;
      tdim = d_image(i) >> IMG2BITS;
      bdim = b.image >> IMG2BITS;
      xdim[2] = IMGMAX + tdim - bdim;

      d_xcmimage(i) = (xdim[2] << IMG2BITS) | (xdim[1] << IMGBITS) | xdim[0];
    }
  );
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  Kokkos::Profiling::pushRegion("rigid/small sort");
  // TODO: check if correct
  if (!setupflag)
    error->all(FLERR, "kk sort before setup");

  auto space = LMPDeviceType();
  auto tmp_bodytag   = Kokkos::create_mirror_view_and_copy(space, d_bodytag);
  auto tmp_bodyown   = Kokkos::create_mirror_view_and_copy(space, d_bodyown);
  auto tmp_xcmimage  = Kokkos::create_mirror_view_and_copy(space, d_xcmimage);
  auto tmp_displace  = Kokkos::create_mirror_view_and_copy(space, d_displace);
  auto tmp_vatom     = Kokkos::create_mirror_view_and_copy(space, d_vatom);
  auto tmp_atom2body = Kokkos::create_mirror_view_and_copy(space, d_atom2body);

  Sorter.sort(space, tmp_bodytag);
  Sorter.sort(space, tmp_bodyown);
  Sorter.sort(space, tmp_xcmimage);
  Sorter.sort(space, tmp_displace);
  Sorter.sort(space, tmp_vatom);
  Sorter.sort(space, tmp_atom2body);

  Kokkos::deep_copy(d_bodytag  , tmp_bodytag  );
  Kokkos::deep_copy(d_bodyown  , tmp_bodyown  );
  Kokkos::deep_copy(d_xcmimage , tmp_xcmimage );
  Kokkos::deep_copy(d_displace , tmp_displace );
  Kokkos::deep_copy(d_vatom    , tmp_vatom    );
  Kokkos::deep_copy(d_atom2body, tmp_atom2body);

  auto d_body = this->d_body;
  auto d_bodyown = this->d_bodyown;
  int nlocal = atom->nlocal;

  Kokkos::parallel_for(
    "fix rigid/small update body.ilocal after sort",
    Range1D(0, nlocal),
    KOKKOS_LAMBDA(const int i){
      if (d_bodyown(i) < 0) return;
      d_body(d_bodyown(i)).ilocal = i;
    }
  );
  Kokkos::Profiling::popRegion();
}


/* ----------------------------------------------------------------------
   zero linear momentum of each rigid body
   set Vcm to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_momentum()
{
  if (!setupflag) {
    FixRigidSmall::zero_momentum();
    return;
  }
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;
  auto nghost_body = this->nghost_body;
  Kokkos::parallel_for(
      "fix rigid/small zero momentum",
      Range1D(0, nlocal_body+nghost_body),
      KOKKOS_LAMBDA(const int ibody){
        double *vcm = d_body(ibody).vcm;
        vcm[0] = vcm[1] = vcm[2] = 0.0;
      }
  );

  // forward communicate vcm to all ghost copies

  commflag = FINAL;
  commKK->forward_comm_device<DeviceType>(this,10);

  // set velocity of atoms in rigid bodues

  evflag = 0;
  set_xv_kokkos(0);
}

/* ----------------------------------------------------------------------
   zero angular momentum of each rigid body
   set angmom/omega to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_rotation()
{
  if (!setupflag) {
    FixRigidSmall::zero_rotation();
    return;
  }
  auto d_body = this->d_body;
  Kokkos::parallel_for(
    "fix rigid/small zero rotation",
    Range1D(0, nlocal_body+nghost_body),
    KOKKOS_LAMBDA(const int ibody){
      double *angmom = d_body(ibody).angmom;
      angmom[0] = angmom[1] = angmom[2] = 0.0;
      double *omega = d_body(ibody).omega;
      omega[0] = omega[1] = omega[2] = 0.0;
    }
  );

  // forward communicate of omega to all ghost copies

  commflag = FINAL;
  commKK->forward_comm_device<DeviceType>(this,10);

  // set velocity of atoms in rigid bodues

  evflag = 0;
  set_xv_kokkos(0);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void *FixRigidSmallKokkos<DeviceType>::extract(const char *str, int &dim)
{
  dim = 0;

  if (strcmp(str,"body") == 0) {
    if (!setupflag) return nullptr;
    dim = 1;
    auto h_atom2body = Kokkos::create_mirror_view_and_copy(LMPHostType(), d_atom2body);
    for(int i = 0; i < atom->nlocal + atom->nghost; i++){
      atom2body[i] = h_atom2body(i);
    }
    return atom2body;
  }

  if (strcmp(str,"onemol") == 0) {
    error->all(FLERR, "onemol not implemented");
    dim = 0;
    return onemols;
  }

  // return vector of rigid body masses, for owned+ghost bodies
  // used by granular pair styles, indexed by atom2body

  if (strcmp(str,"masstotal") == 0) {
    if (!setupflag) return nullptr;
    dim = 1;

    copy_body_host();
    if (nmax_mass < nmax_body) {
      memory->destroy(mass_body);
      nmax_mass = nmax_body;
      memory->create(mass_body,nmax_mass,"rigid:mass_body");
    }

    int n = nlocal_body + nghost_body;
    for (int i = 0; i < n; i++)
      mass_body[i] = body[i].mass;

    return mass_body;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   return translational KE for all rigid bodies
   KE = 1/2 M Vcm^2
   sum local body results across procs
------------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::extract_ke()
{
  if (!setupflag) {
    return FixRigidSmall::extract_ke();
  }

  auto d_body = this->d_body;

  double ke = 0.0;
  Kokkos::parallel_reduce(
    "fix rigid/small ke",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int i, double &ke){
      double *vcm = d_body(i).vcm;
      ke += d_body(i).mass * (vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);
    },
    ke
  );

  double keall;
  MPI_Allreduce(&ke,&keall,1,MPI_DOUBLE,MPI_SUM,world);

  return 0.5*keall;
}

/* ----------------------------------------------------------------------
   return rotational KE for all rigid bodies
   Erotational = 1/2 I wbody^2
------------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::extract_erotational()
{
  if (!setupflag) {
    return FixRigidSmall::extract_erotational();
  }

  double erotate = 0.0;
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;
  Kokkos::parallel_reduce(
    "fix rigid/small erotational",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int i, double &erotate){
      double wbody[3],rot[3][3];
      double *inertia;

      // for Iw^2 rotational term, need wbody = angular velocity in body frame
      // not omega = angular velocity in space frame

      inertia = d_body(i).inertia;
      MathExtraKokkos::quat_to_mat(d_body(i).quat,rot);
      MathExtraKokkos::transpose_matvec(rot,d_body(i).angmom,wbody);
      if (inertia[0] == 0.0) wbody[0] = 0.0;
      else wbody[0] /= inertia[0];
      if (inertia[1] == 0.0) wbody[1] = 0.0;
      else wbody[1] /= inertia[1];
      if (inertia[2] == 0.0) wbody[2] = 0.0;
      else wbody[2] /= inertia[2];

      erotate += inertia[0]*wbody[0]*wbody[0] + inertia[1]*wbody[1]*wbody[1] +
        inertia[2]*wbody[2]*wbody[2];
    },
    erotate
  );

  double erotateall;
  MPI_Allreduce(&erotate,&erotateall,1,MPI_DOUBLE,MPI_SUM,world);

  return 0.5*erotateall;
}

/* ----------------------------------------------------------------------
   return temperature of collection of rigid bodies
   non-active DOF are removed by fflag/tflag and in tfactor
------------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::compute_scalar()
{
  if (!setupflag) {
    return FixRigidSmall::compute_scalar();
  }

  double t = 0.0;
  auto d_body = this->d_body;
  auto nlocal_body = this->nlocal_body;

  Kokkos::parallel_reduce(
    "fix rigid/small compute scalar",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int i, double &t) {
      double wbody[3],rot[3][3];
      double *vcm,*inertia;

      vcm = d_body(i).vcm;
      t += d_body(i).mass * (vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);

      // for Iw^2 rotational term, need wbody = angular velocity in body frame
      // not omega = angular velocity in space frame

      inertia = d_body(i).inertia;
      MathExtraKokkos::quat_to_mat(d_body(i).quat,rot);
      MathExtraKokkos::transpose_matvec(rot,d_body(i).angmom,wbody);
      if (inertia[0] == 0.0) wbody[0] = 0.0;
      else wbody[0] /= inertia[0];
      if (inertia[1] == 0.0) wbody[1] = 0.0;
      else wbody[1] /= inertia[1];
      if (inertia[2] == 0.0) wbody[2] = 0.0;
      else wbody[2] /= inertia[2];

      t += inertia[0]*wbody[0]*wbody[0] + inertia[1]*wbody[1]*wbody[1] +
        inertia[2]*wbody[2]*wbody[2];
    },
    t
  );

  double tall;
  MPI_Allreduce(&t,&tall,1,MPI_DOUBLE,MPI_SUM,world);

  double tfactor = force->mvv2e / ((6.0*nbody - nlinear) * force->boltz);
  tall *= tfactor;
  return tall;
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_body_host(){
  Kokkos::Profiling::pushRegion("rigid/small copy body host");
  auto h_body = Kokkos::create_mirror_view_and_copy(LMPHostType(), d_body);
  for(int ibody = 0; ibody < nlocal_body + nghost_body; ibody++){
    copy_body(&body[ibody], &h_body(ibody));
  }
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_body_device(){
  Kokkos::Profiling::pushRegion("rigid/small copy body device");
  auto h_body = Kokkos::create_mirror_view(d_body);
  for(int ibody = 0; ibody < nlocal_body + nghost_body; ibody++){
    copy_body(&h_body(ibody), &body[ibody]);
  }
  Kokkos::deep_copy(d_body, h_body);
  Kokkos::Profiling::popRegion();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::v_tally(EV_FLOAT &ev, int i, double vtot[6], double r[3], double f[3], double center[3]) const{
  v_tally(ev, i, vtot);

  if (cvflag_atom) {
    const double ri0[3] = {
      r[0]-center[0],
      r[1]-center[1],
      r[2]-center[2],
    };
    cvatom[i][0] += ri0[0]*f[0];
    cvatom[i][1] += ri0[1]*f[1];
    cvatom[i][2] += ri0[2]*f[2];
    cvatom[i][3] += ri0[0]*f[1];
    cvatom[i][4] += ri0[0]*f[2];
    cvatom[i][5] += ri0[1]*f[2];
    cvatom[i][6] += ri0[1]*f[0];
    cvatom[i][7] += ri0[2]*f[0];
    cvatom[i][8] += ri0[2]*f[1];
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixRigidSmallKokkos<DeviceType>::v_tally(EV_FLOAT &ev, int i, double vtot[6]) const{
  if (vflag_global) {
    ev.v[0] += vtot[0];
    ev.v[1] += vtot[1];
    ev.v[2] += vtot[2];
    ev.v[3] += vtot[3];
    ev.v[4] += vtot[4];
    ev.v[5] += vtot[5];
  }

  if (vflag_atom) {
    vatom[i][0] += vtot[0];
    vatom[i][1] += vtot[1];
    vatom[i][2] += vtot[2];
    vatom[i][3] += vtot[3];
    vatom[i][4] += vtot[4];
    vatom[i][5] += vtot[5];
  }
}

namespace LAMMPS_NS {
template class FixRigidSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidSmallKokkos<LMPHostType>;
#endif
}

