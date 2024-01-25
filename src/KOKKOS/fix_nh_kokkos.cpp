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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "fix_nh_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "domain_kokkos.h"
#include "error.h"
#include "fix_deform.h"
#include "force.h"
#include "irregular.h"
#include "kspace.h"
#include "memory_kokkos.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double DELTAFLIP = 0.1;
static constexpr double TILTMAX = 1.5;

enum{NOBIAS,BIAS};
enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

/* ----------------------------------------------------------------------
   NVT,NPH,NPT integrators for improved Nose-Hoover equations of motion
 ---------------------------------------------------------------------- */

template<class DeviceType>
FixNHKokkos<DeviceType>::FixNHKokkos(LAMMPS *lmp, int narg, char **arg) : FixNH(lmp, narg, arg)
{
  kokkosable = 1;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNHKokkos<DeviceType>::init()
{
  FixNH::init();

  atomKK->k_mass.modify<LMPHostType>();
  atomKK->k_mass.sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNHKokkos<DeviceType>::setup(int /*vflag*/)
{
  // tdof needed by compute_temp_target()

  t_current = temperature->compute_scalar();
  tdof = temperature->dof;

  // t_target is needed by NPH and NPT in compute_scalar()
  // If no thermostat or using fix nphug,
  // t_target must be defined by other means.

  if (tstat_flag && strcmp(style,"nphug") != 0) {
    compute_temp_target();
  } else if (pstat_flag) {

    // t0 = reference temperature for masses
    // cannot be done in init() b/c temperature cannot be called there
    // is b/c Modify::init() inits computes after fixes due to dof dependence
    // guesstimate a unit-dependent t0 if actual T = 0.0
    // if it was read in from a restart file, leave it be

    if (t0 == 0.0) {
      atomKK->sync(temperature->execution_space,temperature->datamask_read);
      t0 = temperature->compute_scalar();
      atomKK->modified(temperature->execution_space,temperature->datamask_modify);
      if (t0 == 0.0) {
        if (strcmp(update->unit_style,"lj") == 0) t0 = 1.0;
        else t0 = 300.0;
      }
    }
    t_target = t0;
  }

  if (pstat_flag) compute_press_target();

  atomKK->sync(temperature->execution_space,temperature->datamask_read);
  t_current = temperature->compute_scalar();
  atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  tdof = temperature->dof;

  if (pstat_flag) {
    //atomKK->sync(pressure->execution_space,pressure->datamask_read);
    //atomKK->modified(pressure->execution_space,pressure->datamask_modify);
    if (pstyle == ISO) pressure->compute_scalar();
    else pressure->compute_vector();
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  // masses and initial forces on thermostat variables

  if (tstat_flag) {
    eta_mass[0] = tdof * boltz * t_target / (t_freq*t_freq);
    for (int ich = 1; ich < mtchain; ich++)
      eta_mass[ich] = boltz * t_target / (t_freq*t_freq);
    for (int ich = 1; ich < mtchain; ich++) {
      eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1] -
                         boltz * t_target) / eta_mass[ich];
    }
  }

  // masses and initial forces on barostat variables

  if (pstat_flag) {
    double kt = boltz * t_target;
    double nkt = (atom->natoms + 1) * kt;

    for (int i = 0; i < 3; i++)
      if (p_flag[i])
        omega_mass[i] = nkt/(p_freq[i]*p_freq[i]);

    if (pstyle == TRICLINIC) {
      for (int i = 3; i < 6; i++)
        if (p_flag[i]) omega_mass[i] = nkt/(p_freq[i]*p_freq[i]);
    }

  // masses and initial forces on barostat thermostat variables

    if (mpchain) {
      etap_mass[0] = boltz * t_target / (p_freq_max*p_freq_max);
      for (int ich = 1; ich < mpchain; ich++)
        etap_mass[ich] = boltz * t_target / (p_freq_max*p_freq_max);
      for (int ich = 1; ich < mpchain; ich++)
        etap_dotdot[ich] =
          (etap_mass[ich-1]*etap_dot[ich-1]*etap_dot[ich-1] -
           boltz * t_target) / etap_mass[ich];
    }

  }
}

/* ----------------------------------------------------------------------
   1st half of Verlet update
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNHKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  // update eta_press_dot

  if (pstat_flag && mpchain) nhc_press_integrate();

  // update eta_dot

  if (tstat_flag) {
    compute_temp_target();
    nhc_temp_integrate();
  }

  // need to recompute pressure to account for change in KE
  // t_current is up-to-date, but compute_temperature is not
  // compute appropriately coupled elements of mvv_current

  if (pstat_flag) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
    //atomKK->sync(pressure->execution_space,pressure->datamask_read);
    //atomKK->modified(pressure->execution_space,pressure->datamask_modify);
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  if (pstat_flag) {
    compute_press_target();
    nh_omega_dot();
    nh_v_press();
  }

  nve_v();

  // remap simulation box by 1/2 step

  if (pstat_flag) remap();

  nve_x();

  // remap simulation box by 1/2 step
  // redo KSpace coeffs since volume has changed

  if (pstat_flag) {
    remap();
    if (kspace_flag) force->kspace->setup();
  }
}

/* ----------------------------------------------------------------------
   2nd half of Verlet update
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNHKokkos<DeviceType>::final_integrate()
{
  nve_v();

  // re-compute temp before nh_v_press()
  // only needed for temperature computes with BIAS on reneighboring steps:
  //   b/c some biases store per-atom values (e.g. temp/profile)
  //   per-atom values are invalid if reneigh/comm occurred
  //     since temp->compute() in initial_integrate()

  if (which == BIAS && neighbor->ago == 0) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    t_current = temperature->compute_scalar();
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  }

  if (pstat_flag) nh_v_press();

  // compute new T,P after velocities rescaled by nh_v_press()
  // compute appropriately coupled elements of mvv_current

  atomKK->sync(temperature->execution_space,temperature->datamask_read);
  t_current = temperature->compute_scalar();
  atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  tdof = temperature->dof;

  if (pstat_flag) {
    //atomKK->sync(pressure->execution_space,pressure->datamask_read);
    //atomKK->modified(pressure->execution_space,pressure->datamask_modify);
    if (pstyle == ISO) pressure->compute_scalar();
    else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  if (pstat_flag) nh_omega_dot();

  // update eta_dot
  // update eta_press_dot

  if (tstat_flag) nhc_temp_integrate();
  if (pstat_flag && mpchain) nhc_press_integrate();
}

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or dilate group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNHKokkos<DeviceType>::remap()
{
  double oldlo,oldhi;
  double expfac;

  int nlocal = atom->nlocal;
  double *h = domain->h;

  // omega is not used, except for book-keeping

  for (int i = 0; i < 6; i++) omega[i] += dto*omega_dot[i];

  // convert pertinent atoms and rigid bodies to lamda coords

  domainKK->x2lamda(nlocal);
  //if (allremap) domainKK->x2lamda(nlocal);
  //else {
  //  for (i = 0; i < nlocal; i++)
  //    if (mask[i] & dilate_group_bit)
  //      domain->x2lamda(x[i],x[i]);
  //}

  if (rfix.size() > 0)
    error->all(FLERR,"Cannot (yet) use rigid bodies with fix nh and Kokkos");
  // for (auto &ifix : rfix) ifix->deform(0);

  // reset global and local box to new size/shape

  // this operation corresponds to applying the
  // translate and scale operations
  // corresponding to the solution of the following ODE:
  //
  // h_dot = omega_dot * h
  //
  // where h_dot, omega_dot and h are all upper-triangular
  // 3x3 tensors. In Voigt notation, the elements of the
  // RHS product tensor are:
  // h_dot = [0*0, 1*1, 2*2, 1*3+3*2, 0*4+5*3+4*2, 0*5+5*1]
  //
  // Ordering of operations preserves time symmetry.

  double dto2 = dto/2.0;
  double dto4 = dto/4.0;
  double dto8 = dto/8.0;

  // off-diagonal components, first half

  if (pstyle == TRICLINIC) {

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }

    if (p_flag[3]) {
      expfac = exp(dto4*omega_dot[1]);
      h[3] *= expfac;
      h[3] += dto2*(omega_dot[3]*h[2]);
      h[3] *= expfac;
    }

    if (p_flag[5]) {
      expfac = exp(dto4*omega_dot[0]);
      h[5] *= expfac;
      h[5] += dto2*(omega_dot[5]*h[1]);
      h[5] *= expfac;
    }

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }
  }

  // scale diagonal components
  // scale tilt factors with cell, if set

  if (p_flag[0]) {
    oldlo = domain->boxlo[0];
    oldhi = domain->boxhi[0];
    expfac = exp(dto*omega_dot[0]);
    domain->boxlo[0] = (oldlo-fixedpoint[0])*expfac + fixedpoint[0];
    domain->boxhi[0] = (oldhi-fixedpoint[0])*expfac + fixedpoint[0];
  }

  if (p_flag[1]) {
    oldlo = domain->boxlo[1];
    oldhi = domain->boxhi[1];
    expfac = exp(dto*omega_dot[1]);
    domain->boxlo[1] = (oldlo-fixedpoint[1])*expfac + fixedpoint[1];
    domain->boxhi[1] = (oldhi-fixedpoint[1])*expfac + fixedpoint[1];
    if (scalexy) h[5] *= expfac;
  }

  if (p_flag[2]) {
    oldlo = domain->boxlo[2];
    oldhi = domain->boxhi[2];
    expfac = exp(dto*omega_dot[2]);
    domain->boxlo[2] = (oldlo-fixedpoint[2])*expfac + fixedpoint[2];
    domain->boxhi[2] = (oldhi-fixedpoint[2])*expfac + fixedpoint[2];
    if (scalexz) h[4] *= expfac;
    if (scaleyz) h[3] *= expfac;
  }

  // off-diagonal components, second half

  if (pstyle == TRICLINIC) {

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }

    if (p_flag[3]) {
      expfac = exp(dto4*omega_dot[1]);
      h[3] *= expfac;
      h[3] += dto2*(omega_dot[3]*h[2]);
      h[3] *= expfac;
    }

    if (p_flag[5]) {
      expfac = exp(dto4*omega_dot[0]);
      h[5] *= expfac;
      h[5] += dto2*(omega_dot[5]*h[1]);
      h[5] *= expfac;
    }

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }

  }

  domain->yz = h[3];
  domain->xz = h[4];
  domain->xy = h[5];

  // tilt factor to cell length ratio can not exceed TILTMAX in one step

  if (domain->yz < -TILTMAX*domain->yprd ||
      domain->yz > TILTMAX*domain->yprd ||
      domain->xz < -TILTMAX*domain->xprd ||
      domain->xz > TILTMAX*domain->xprd ||
      domain->xy < -TILTMAX*domain->xprd ||
      domain->xy > TILTMAX*domain->xprd)
    error->all(FLERR,"Fix npt/nph has tilted box too far in one step - "
               "periodic cell is too far from equilibrium state");

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  domainKK->lamda2x(nlocal);
  //if (allremap) domainKK->lamda2x(nlocal);
  //else {
  //  for (i = 0; i < nlocal; i++)
  //    if (mask[i] & dilate_group_bit)
  //      domain->lamda2x(x[i],x[i]);
  //}

  // for (auto &ifix : rfix) ifix->deform(1);
}

/* ----------------------------------------------------------------------
   perform half-step barostat scaling of velocities
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixNHKokkos<DeviceType>::nh_v_press()
{
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  factor[0] = exp(-dt4*(omega_dot[0]+mtk_term2));
  factor[1] = exp(-dt4*(omega_dot[1]+mtk_term2));
  factor[2] = exp(-dt4*(omega_dot[2]+mtk_term2));

  if (which == BIAS) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    temperature->remove_bias_all();
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  }

  atomKK->sync(execution_space,V_MASK | MASK_MASK);

  copymode = 1;
  if (pstyle == TRICLINIC)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNH_nh_v_press<1> >(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNH_nh_v_press<0> >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);

  if (which == BIAS) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    temperature->restore_bias_all();
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  }

}

template<class DeviceType>
template<int TRICLINIC_FLAG>
KOKKOS_INLINE_FUNCTION
void FixNHKokkos<DeviceType>::operator()(TagFixNH_nh_v_press<TRICLINIC_FLAG>, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,0) *= factor[0];
    v(i,1) *= factor[1];
    v(i,2) *= factor[2];
    if (TRICLINIC_FLAG) {
      v(i,0) += -dthalf*(v(i,1)*omega_dot[5] + v(i,2)*omega_dot[4]);
      v(i,1) += -dthalf*v(i,2)*omega_dot[3];
    }
    v(i,0) *= factor[0];
    v(i,1) *= factor[1];
    v(i,2) *= factor[2];
  }
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixNHKokkos<DeviceType>::nve_v()
{
  atomKK->sync(execution_space,X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);

  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  copymode = 1;
  if (rmass.data())
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNH_nve_v<1> >(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNH_nve_v<0> >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
template<int RMASS>
KOKKOS_INLINE_FUNCTION
void FixNHKokkos<DeviceType>::operator()(TagFixNH_nve_v<RMASS>, const int &i) const {
  if (RMASS) {
    if (mask[i] & groupbit) {
      const F_FLOAT dtfm = dtf / rmass[i];
      v(i,0) += dtfm*f(i,0);
      v(i,1) += dtfm*f(i,1);
      v(i,2) += dtfm*f(i,2);
    }
  } else {
    if (mask[i] & groupbit) {
      const F_FLOAT dtfm = dtf / mass[type[i]];
      v(i,0) += dtfm*f(i,0);
      v(i,1) += dtfm*f(i,1);
      v(i,2) += dtfm*f(i,2);
    }
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixNHKokkos<DeviceType>::nve_x()
{
  atomKK->sync(execution_space,X_MASK | V_MASK | MASK_MASK);
  atomKK->modified(execution_space,X_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  // x update by full step only for atoms in group

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNH_nve_x>(0,nlocal),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNHKokkos<DeviceType>::operator()(TagFixNH_nve_x, const int &i) const {
  if (mask[i] & groupbit) {
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

/* ----------------------------------------------------------------------
   perform half-step thermostat scaling of velocities
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixNHKokkos<DeviceType>::nh_v_temp()
{
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  if (which == BIAS) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    temperature->remove_bias_all();
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  }

  atomKK->sync(execution_space,V_MASK | MASK_MASK);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNH_nh_v_temp>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);

  if (which == BIAS) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    temperature->restore_bias_all();
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNHKokkos<DeviceType>::operator()(TagFixNH_nh_v_temp, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,0) *= factor_eta;
    v(i,1) *= factor_eta;
    v(i,2) *= factor_eta;
  }
}

/* ----------------------------------------------------------------------
  if any tilt ratios exceed limits, set flip = 1 and compute new tilt values
  do not flip in x or y if non-periodic (can tilt but not flip)
    this is b/c the box length would be changed (dramatically) by flip
  if yz tilt exceeded, adjust C vector by one B vector
  if xz tilt exceeded, adjust C vector by one A vector
  if xy tilt exceeded, adjust B vector by one A vector
  check yz first since it may change xz, then xz check comes after
  if any flip occurs, create new box in domain
  image_flip() adjusts image flags due to box shape change induced by flip
  remap() puts atoms outside the new box back into the new box
  perform irregular on atoms in lamda coords to migrate atoms to new procs
  important that image_flip comes before remap, since remap may change
    image flags to new values, making eqs in doc of Domain:image_flip incorrect
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNHKokkos<DeviceType>::pre_exchange()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;

  // flip is only triggered when tilt exceeds 0.5 by DELTAFLIP
  // this avoids immediate re-flipping due to tilt oscillations

  double xtiltmax = (0.5+DELTAFLIP)*xprd;
  double ytiltmax = (0.5+DELTAFLIP)*yprd;

  int flipxy,flipxz,flipyz;
  flipxy = flipxz = flipyz = 0;

  if (domain->yperiodic) {
    if (domain->yz < -ytiltmax) {
      domain->yz += yprd;
      domain->xz += domain->xy;
      flipyz = 1;
    } else if (domain->yz >= ytiltmax) {
      domain->yz -= yprd;
      domain->xz -= domain->xy;
      flipyz = -1;
    }
  }

  if (domain->xperiodic) {
    if (domain->xz < -xtiltmax) {
      domain->xz += xprd;
      flipxz = 1;
    } else if (domain->xz >= xtiltmax) {
      domain->xz -= xprd;
      flipxz = -1;
    }
    if (domain->xy < -xtiltmax) {
      domain->xy += xprd;
      flipxy = 1;
    } else if (domain->xy >= xtiltmax) {
      domain->xy -= xprd;
      flipxy = -1;
    }
  }

  int flip = 0;
  if (flipxy || flipxz || flipyz) flip = 1;

  if (flip) {
    domain->set_global_box();
    domain->set_local_box();

    domainKK->image_flip(flipxy,flipxz,flipyz);

    domainKK->remap_all();

    domainKK->x2lamda(atom->nlocal);
    atomKK->sync(Host,ALL_MASK);
    irregular->migrate_atoms();
    atomKK->modified(Host,ALL_MASK);
    domainKK->lamda2x(atom->nlocal);
  }
}

namespace LAMMPS_NS {
template class FixNHKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNHKokkos<LMPHostType>;
#endif
}
