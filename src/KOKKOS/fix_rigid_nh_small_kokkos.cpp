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
   Contributing author: Trung Dac Nguyen (ORNL) for the host fix rigid/nh/small
   references: Kamberaj et al., J. Chem. Phys. 122, 224114 (2005)
               Miller et al., J Chem Phys. 116, 8649-8659 (2002)
------------------------------------------------------------------------- */

#include "fix_rigid_nh_small_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_extra.h"
#include "math_extra_kokkos.h"
#include "memory.h"
#include "modify.h"
#include "rigid_const.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace RigidConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNHSmallKokkos<DeviceType>::FixRigidNHSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
    FixRigidSmallKokkos<DeviceType>(lmp, narg, arg), w(nullptr), wdti1(nullptr), wdti2(nullptr),
    wdti4(nullptr), q_t(nullptr), q_r(nullptr), eta_t(nullptr), eta_r(nullptr), eta_dot_t(nullptr),
    eta_dot_r(nullptr), f_eta_t(nullptr), f_eta_r(nullptr), q_b(nullptr), eta_b(nullptr),
    eta_dot_b(nullptr), f_eta_b(nullptr), id_temp(nullptr), id_press(nullptr),
    temperature(nullptr), pressure(nullptr)
{
  if (tstat_flag || pstat_flag) this->ecouple_flag = 1;

  // error checks (mirror FixRigidNHSmall)

  if (domain->dimension == 2 && p_flag[2])
    error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);
  if (domain->dimension == 2 && (pcouple == YZ || pcouple == XZ))
    error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);

  if (pcouple == XYZ && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == XYZ && domain->dimension == 3 && p_flag[2] == 0)
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);

  if (p_flag[0] && domain->xperiodic == 0)
    error->all(FLERR, "Cannot use fix {} on a non-periodic dimension", style);
  if (p_flag[1] && domain->yperiodic == 0)
    error->all(FLERR, "Cannot use fix {} on a non-periodic dimension", style);
  if (p_flag[2] && domain->zperiodic == 0)
    error->all(FLERR, "Cannot use fix {} on a non-periodic dimension", style);

  if (pcouple == XYZ && domain->dimension == 3 &&
      (p_start[0] != p_start[1] || p_start[0] != p_start[2] ||
       p_stop[0] != p_stop[1] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[1] || p_period[0] != p_period[2]))
    error->all(FLERR, "Invalid fix {} command pressure settings", style);
  if (pcouple == XYZ && domain->dimension == 2 &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR, "Invalid fix {} command pressure settings", style);
  if (pcouple == XY &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR, "Invalid fix {} command pressure settings", style);
  if (pcouple == YZ &&
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2] ||
       p_period[1] != p_period[2]))
    error->all(FLERR, "Invalid fix {} command pressure settings", style);
  if (pcouple == XZ &&
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[2]))
    error->all(FLERR, "Invalid fix {} command pressure settings", style);

  if (p_flag[0]) this->box_change |= Fix::BOX_CHANGE_X;
  if (p_flag[1]) this->box_change |= Fix::BOX_CHANGE_Y;
  if (p_flag[2]) this->box_change |= Fix::BOX_CHANGE_Z;

  if ((tstat_flag && t_period <= 0.0) ||
      (p_flag[0] && p_period[0] <= 0.0) ||
      (p_flag[1] && p_period[1] <= 0.0) ||
      (p_flag[2] && p_period[2] <= 0.0))
    error->all(FLERR,"Fix {} damping parameters must be > 0.0", style);

  // memory allocation and initialization

  if (tstat_flag || pstat_flag) {
    allocate_chain();
    allocate_order();
  }

  if (tstat_flag) {
    eta_t[0] = eta_r[0] = 0.0;
    eta_dot_t[0] = eta_dot_r[0] = 0.0;
    f_eta_t[0] = f_eta_r[0] = 0.0;

    for (int i = 1; i < t_chain; i++) {
      eta_t[i] = eta_r[i] = 0.0;
      eta_dot_t[i] = eta_dot_r[i] = 0.0;
    }
  }

  if (pstat_flag) {
    epsilon_dot[0] = epsilon_dot[1] = epsilon_dot[2] = 0.0;
    eta_b[0] = eta_dot_b[0] = f_eta_b[0] = 0.0;
    for (int i = 1; i < p_chain; i++)
      eta_b[i] = eta_dot_b[i] = 0.0;
  }

  vol0 = 0.0;
  t0 = 1.0;

  tcomputeflag = 0;
  pcomputeflag = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNHSmallKokkos<DeviceType>::~FixRigidNHSmallKokkos()
{
  if (this->copymode) return;

  if (tstat_flag || pstat_flag) {
    deallocate_chain();
    deallocate_order();
  }

  if (tcomputeflag) modify->delete_compute(id_temp);
  delete[] id_temp;

  if (pstat_flag) {
    if (pcomputeflag) modify->delete_compute(id_press);
    delete[] id_press;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidNHSmallKokkos<DeviceType>::setmask()
{
  // inherit the kk base mask (adds PRE_EXCHANGE for host atom migration)
  return FixRigidSmallKokkos<DeviceType>::setmask();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::init()
{
  FixRigidSmallKokkos<DeviceType>::init();

  // recheck that dilate group has not been deleted

  if (allremap == 0) {
    int idilate = group->find(id_dilate);
    if (idilate == -1) error->all(FLERR,"Fix {} dilate group ID does not exist", style);
    dilate_group_bit = group->bitmask[idilate];
  }

  // initialize thermostats / set timesteps, constants, Yoshida-Suzuki params

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;

  boltz = force->boltz;
  nktv2p = force->nktv2p;
  mvv2e = force->mvv2e;

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // see Table 1 in Kamberaj et al

  if (tstat_flag || pstat_flag) {
    if (t_order == 3) {
      w[0] = 1.0 / (2.0 - pow(2.0, 1.0/3.0));
      w[1] = 1.0 - 2.0*w[0];
      w[2] = w[0];
    } else if (t_order == 5) {
      w[0] = 1.0 / (4.0 - pow(4.0, 1.0/3.0));
      w[1] = w[0];
      w[2] = 1.0 - 4.0 * w[0];
      w[3] = w[0];
      w[4] = w[0];
    }
  }

  if (tcomputeflag) {
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature) {
      error->all(FLERR,"Temperature compute ID {} for fix {} does not exist", id_temp, style);
    } else {
      if (temperature->tempflag == 0)
        error->all(FLERR, "Compute ID {} for fix {} does not compute a temperature", id_temp, style);
    }
  }

  if (pstat_flag) {
    if (domain->triclinic)
      error->all(FLERR,"Fix {} does not yet allow triclinic box", style);

    // ensure no conflict with fix deform

    for (const auto &ifix : modify->get_fix_by_style("^deform")) {
      auto *deform = dynamic_cast<FixDeform *>(ifix);
      if (deform) {
        int *dimflag = deform->dimflag;
        if ((p_flag[0] && dimflag[0]) || (p_flag[1] && dimflag[1]) ||
            (p_flag[2] && dimflag[2]))
          error->all(FLERR,"Cannot use fix {} and fix deform on same component of stress tensor",
                     style);
      }
    }

    // set frequency

    p_freq_max = 0.0;
    p_freq_max = MAX(p_freq[0],p_freq[1]);
    p_freq_max = MAX(p_freq_max,p_freq[2]);

    // tally the number of dimensions that are barostatted
    // set initial volume and reference cell, if not already done

    pdim = p_flag[0] + p_flag[1] + p_flag[2];
    if (vol0 == 0.0) {
      if (domain->dimension == 2) vol0 = domain->xprd * domain->yprd;
      else vol0 = domain->xprd * domain->yprd * domain->zprd;
    }

    // set pressure compute ptr

    pressure = modify->get_compute_by_id(id_press);
    if (!pressure) {
      error->all(FLERR,"Pressure compute ID {} for fix {} does not exist", id_press, style);
    } else {
      if (pressure->pressflag == 0)
        error->all(FLERR,"Compute ID {} for fix {} does not compute pressure", id_press, style);
    }

    // detect if any rigid fixes exist so rigid bodies move on remap

    rfix.clear();
    for (const auto &ifix : modify->get_fix_list())
      if (ifix->rigid_flag) rfix.push_back(ifix);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::setup(int vflag)
{
  // host setup populates body[] and the per-atom arrays
  atomKK->sync(Host, datamask_read);
  FixRigidSmall::setup(vflag);
  atomKK->modified(Host, datamask_modify);
  // reconcile the legacy-host atom:x/v writes into the fix's execution_space so
  // the later ModifyKokkos host-modify does not trip the TransformView guard in
  // a mixed/single build (no-op in a full-double build); see the base setup()
  atomKK->sync(execution_space, datamask_modify);

  // --- Nose-Hoover scalar setup (host; reads/writes host body[]) ---

  compute_dof();

  double mbody[3];
  akin_t = akin_r = 0.0;
  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];
    MathExtra::transpose_matvec(b->ex_space,b->ey_space,b->ez_space,b->angmom,mbody);
    MathExtra::quatvec(b->quat,mbody,b->conjqm);
    b->conjqm[0] *= 2.0;
    b->conjqm[1] *= 2.0;
    b->conjqm[2] *= 2.0;
    b->conjqm[3] *= 2.0;

    if (tstat_flag || pstat_flag) {
      akin_t += b->mass*(b->vcm[0]*b->vcm[0] + b->vcm[1]*b->vcm[1] + b->vcm[2]*b->vcm[2]);
      akin_r += b->angmom[0]*b->omega[0] + b->angmom[1]*b->omega[1] + b->angmom[2]*b->omega[2];
    }
  }

  if (tstat_flag || pstat_flag) {
    double ke[2],keall[2];
    ke[0] = akin_t;
    ke[1] = akin_r;
    MPI_Allreduce(ke,keall,2,MPI_DOUBLE,MPI_SUM,world);
    akin_t = keall[0];
    akin_r = keall[1];
  }

  if (tstat_flag) compute_temp_target();
  else if (pstat_flag) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    t0 = temperature->compute_scalar();
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
    atomKK->sync(execution_space,temperature->datamask_modify);
    if (t0 == 0.0) {
      if (strcmp(update->unit_style,"lj") == 0) t0 = 1.0;
      else t0 = 300.0;
    }
    t_target = t0;
  }

  if (pstat_flag) {
    compute_press_target();

    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    atomKK->sync(pressure->execution_space,pressure->datamask_read);
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
    atomKK->modified(pressure->execution_space,pressure->datamask_modify);
    atomKK->sync(execution_space,temperature->datamask_modify);
    atomKK->sync(execution_space,pressure->datamask_modify);

    couple();
    pressure->addstep(update->ntimestep+1);
  }

  // initialize thermostat/barostat settings

  double kt, t_mass, tb_mass;
  kt = boltz * t_target;

  if (tstat_flag) {
    t_mass = kt / (t_freq*t_freq);
    q_t[0] = nf_t * t_mass;
    q_r[0] = nf_r * t_mass;
    for (int i = 1; i < t_chain; i++)
      q_t[i] = q_r[i] = t_mass;

    for (int i = 1; i < t_chain; i++) {
      f_eta_t[i] = (q_t[i-1] * eta_dot_t[i-1] * eta_dot_t[i-1] - kt)/q_t[i];
      f_eta_r[i] = (q_r[i-1] * eta_dot_r[i-1] * eta_dot_r[i-1] - kt)/q_r[i];
    }
  }

  int dimension = domain->dimension;

  if (pstat_flag) {
    for (int i = 0; i < 3; i++)
      if (p_flag[i]) {
        epsilon_mass[i] = (g_f + dimension) * kt / (p_freq[i]*p_freq[i]);
        epsilon[i] = log(vol0)/dimension;
      }

    tb_mass = kt / (p_freq_max * p_freq_max);
    q_b[0] = dimension * dimension * tb_mass;
    for (int i = 1; i < p_chain; i++) {
      q_b[i] = tb_mass;
      f_eta_b[i] = (q_b[i] * eta_dot_b[i-1] * eta_dot_b[i-1] - kt)/q_b[i];
    }
  }

  if (tstat_flag || pstat_flag) {
    for (int i = 0; i < t_order; i++) {
      wdti1[i] = w[i] * dtv / t_iter;
      wdti2[i] = wdti1[i] / 2.0;
      wdti4[i] = wdti1[i] / 4.0;
    }
  }

  if (pstat_flag) {
    compute_press_target();
    nh_epsilon_dot();
  }

  // --- push host state to device and enable device comm/sort for the run ---

  setup_device_push();
}

/* ----------------------------------------------------------------------
   perform preforce velocity Verlet integration
   see Kamberaj paper for step references
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::initial_integrate(int vflag)
{
  double tmp,scale_r,scale_t[3],scale_v[3];

  scale_t[0] = scale_t[1] = scale_t[2] = 1.0;
  scale_v[0] = scale_v[1] = scale_v[2] = 1.0;
  scale_r = 1.0;

  if (tstat_flag) {
    tmp = exp(-dtq * eta_dot_t[0]);
    scale_t[0] = scale_t[1] = scale_t[2] = tmp;
    tmp = exp(-dtq * eta_dot_r[0]);
    scale_r = tmp;
  }

  if (pstat_flag) {
    scale_t[0] *= exp(-dtq * (epsilon_dot[0] + mtk_term2));
    scale_t[1] *= exp(-dtq * (epsilon_dot[1] + mtk_term2));
    scale_t[2] *= exp(-dtq * (epsilon_dot[2] + mtk_term2));
    scale_r *= exp(-dtq * (pdim * mtk_term2));

    tmp = dtq * epsilon_dot[0];
    scale_v[0] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[1];
    scale_v[1] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[2];
    scale_v[2] = dtv * exp(tmp) * maclaurin_series(tmp);
  }

  // update xcm, vcm, quat, conjqm and angmom of bodies on the device

  nh_initial_integrate_bodies(scale_r, scale_t, scale_v);

  // forward communicate updated info of all bodies

  commflag = INITIAL;
  commKK->template forward_comm_device<DeviceType>(this,29);

  // accumulate translational and rotational kinetic energies

  if (tstat_flag || pstat_flag) nh_akin(akin_t, akin_r);

  // compute target temperature, update thermostat chains

  if (tstat_flag) {
    compute_temp_target();
    if (dynamic) { copy_body_host(); compute_dof(); }
    nhc_temp_integrate();
  }

  if (pstat_flag) nhc_press_integrate();

  // virial setup before call to set_xv

  this->v_init(vflag);

  // remap simulation box by 1/2 step

  if (pstat_flag) remap();

  // set coords/orient and velocity/rotation of atoms in rigid bodies

  set_xv_kokkos(1);

  // remap simulation box by full step; redo KSpace coeffs since volume changed

  if (pstat_flag) {
    remap();
    if (kspace_flag) force->kspace->setup();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::final_integrate()
{
  double tmp,scale_t[3],scale_r;

  scale_t[0] = scale_t[1] = scale_t[2] = 1.0;
  scale_r = 1.0;

  if (tstat_flag) {
    tmp = exp(-1.0 * dtq * eta_dot_t[0]);
    scale_t[0] = scale_t[1] = scale_t[2] = tmp;
    scale_r = exp(-1.0 * dtq * eta_dot_r[0]);
  }

  if (pstat_flag) {
    scale_t[0] *= exp(-dtq * (epsilon_dot[0] + mtk_term2));
    scale_t[1] *= exp(-dtq * (epsilon_dot[1] + mtk_term2));
    scale_t[2] *= exp(-dtq * (epsilon_dot[2] + mtk_term2));
    scale_r *= exp(-dtq * (pdim * mtk_term2));
  }

  // late calculation of forces and torques (if requested)

  if (!earlyflag) compute_forces_and_torques_kokkos();

  // update vcm and angmom, recompute omega, on the device

  nh_final_integrate_bodies(scale_r, scale_t);

  // forward communicate updated info of all bodies

  commflag = FINAL;
  commKK->template forward_comm_device<DeviceType>(this,10);

  // accumulate translational and rotational kinetic energies

  if (pstat_flag) nh_akin(akin_t, akin_r);

  // set velocity/rotation of atoms in rigid bodies (virial set in initial_integrate)

  set_xv_kokkos(0);

  // compute current temperature

  if (tcomputeflag) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    t_current = temperature->compute_scalar();
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
    atomKK->sync(execution_space,temperature->datamask_modify);
  }

  // compute current and target pressures, update epsilon_dot

  if (pstat_flag) {
    atomKK->sync(temperature->execution_space,temperature->datamask_read);
    atomKK->sync(pressure->execution_space,pressure->datamask_read);
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    atomKK->modified(temperature->execution_space,temperature->datamask_modify);
    atomKK->modified(pressure->execution_space,pressure->datamask_modify);
    atomKK->sync(execution_space,temperature->datamask_modify);
    atomKK->sync(execution_space,pressure->datamask_modify);

    couple();
    pressure->addstep(update->ntimestep+1);

    compute_press_target();

    nh_epsilon_dot();
  }
}

/* ----------------------------------------------------------------------
   per-body initial Nose-Hoover update on the device
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::nh_initial_integrate_bodies(
    const double scale_r, const double scale_t_in[3], const double scale_v_in[3])
{
  auto d_body = this->d_body;
  const double l_dtf = dtf;
  const double l_dtv = dtv;
  const double l_dtq = dtq;
  const double dtf2 = dtf * 2.0;
  const int scaleflag = (tstat_flag || pstat_flag) ? 1 : 0;
  const int pflag = pstat_flag ? 1 : 0;
  const double scale_r_l = scale_r;
  double scale_t[3] = {scale_t_in[0], scale_t_in[1], scale_t_in[2]};
  double scale_v[3] = {scale_v_in[0], scale_v_in[1], scale_v_in[2]};

  Kokkos::parallel_for("fix rigid/nh/small initial body",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = d_body(ibody);

      // step 1.1 - update vcm by 1/2 step

      double dtfm = l_dtf / b.mass;
      b.vcm[0] += dtfm * b.fcm[0];
      b.vcm[1] += dtfm * b.fcm[1];
      b.vcm[2] += dtfm * b.fcm[2];

      if (scaleflag) {
        b.vcm[0] *= scale_t[0];
        b.vcm[1] *= scale_t[1];
        b.vcm[2] *= scale_t[2];
      }

      // step 1.2 - update xcm by full step

      if (!pflag) {
        b.xcm[0] += l_dtv * b.vcm[0];
        b.xcm[1] += l_dtv * b.vcm[1];
        b.xcm[2] += l_dtv * b.vcm[2];
      } else {
        b.xcm[0] += scale_v[0] * b.vcm[0];
        b.xcm[1] += scale_v[1] * b.vcm[1];
        b.xcm[2] += scale_v[2] * b.vcm[2];
      }

      // step 1.3 - apply torque (body coords) to quaternion momentum
      // rotational quantities are computed in KK_FLOAT; the orientation and
      // momentum quaternions (quat, conjqm, fquat) stay double

      KK_FLOAT ex[3] = {(KK_FLOAT)b.ex_space[0],(KK_FLOAT)b.ex_space[1],(KK_FLOAT)b.ex_space[2]};
      KK_FLOAT ey[3] = {(KK_FLOAT)b.ey_space[0],(KK_FLOAT)b.ey_space[1],(KK_FLOAT)b.ey_space[2]};
      KK_FLOAT ez[3] = {(KK_FLOAT)b.ez_space[0],(KK_FLOAT)b.ez_space[1],(KK_FLOAT)b.ez_space[2]};
      KK_FLOAT torque[3] = {(KK_FLOAT)b.torque[0],(KK_FLOAT)b.torque[1],(KK_FLOAT)b.torque[2]};
      KK_FLOAT inertia[3] = {(KK_FLOAT)b.inertia[0],(KK_FLOAT)b.inertia[1],(KK_FLOAT)b.inertia[2]};
      KK_FLOAT tbody[3], mbody[3], angmom[3], omega[3];
      double fquat[4];
      MathExtraKokkos::transpose_matvec(ex,ey,ez,torque,tbody);
      MathExtraKokkos::quatvec(b.quat,tbody,fquat);

      b.conjqm[0] += dtf2 * fquat[0];
      b.conjqm[1] += dtf2 * fquat[1];
      b.conjqm[2] += dtf2 * fquat[2];
      b.conjqm[3] += dtf2 * fquat[3];

      if (scaleflag) {
        b.conjqm[0] *= scale_r_l;
        b.conjqm[1] *= scale_r_l;
        b.conjqm[2] *= scale_r_l;
        b.conjqm[3] *= scale_r_l;
      }

      // step 1.4 to 1.13 - use no_squish rotate to update p and q

      MathExtraKokkos::no_squish_rotate(3,b.conjqm,b.quat,inertia,l_dtq);
      MathExtraKokkos::no_squish_rotate(2,b.conjqm,b.quat,inertia,l_dtq);
      MathExtraKokkos::no_squish_rotate(1,b.conjqm,b.quat,inertia,l_dtv);
      MathExtraKokkos::no_squish_rotate(2,b.conjqm,b.quat,inertia,l_dtq);
      MathExtraKokkos::no_squish_rotate(3,b.conjqm,b.quat,inertia,l_dtq);

      // update exyz_space, transform p back to angmom, update angular velocity

      MathExtraKokkos::q_to_exyz(b.quat,ex,ey,ez);
      MathExtraKokkos::invquatvec(b.quat,b.conjqm,mbody);
      MathExtraKokkos::matvec(ex,ey,ez,mbody,angmom);

      angmom[0] *= 0.5;
      angmom[1] *= 0.5;
      angmom[2] *= 0.5;

      MathExtraKokkos::angmom_to_omega(angmom,ex,ey,ez,inertia,omega);

      // store the updated rotational state back to the (double) body
      b.ex_space[0]=ex[0]; b.ex_space[1]=ex[1]; b.ex_space[2]=ex[2];
      b.ey_space[0]=ey[0]; b.ey_space[1]=ey[1]; b.ey_space[2]=ey[2];
      b.ez_space[0]=ez[0]; b.ez_space[1]=ez[1]; b.ez_space[2]=ez[2];
      b.angmom[0]=angmom[0]; b.angmom[1]=angmom[1]; b.angmom[2]=angmom[2];
      b.omega[0]=omega[0]; b.omega[1]=omega[1]; b.omega[2]=omega[2];
    });
}

/* ----------------------------------------------------------------------
   per-body final Nose-Hoover update on the device
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::nh_final_integrate_bodies(
    const double scale_r, const double scale_t_in[3])
{
  auto d_body = this->d_body;
  const double l_dtf = dtf;
  const double dtf2 = dtf * 2.0;
  const int scaleflag = (tstat_flag || pstat_flag) ? 1 : 0;
  const double scale_r_l = scale_r;
  double scale_t[3] = {scale_t_in[0], scale_t_in[1], scale_t_in[2]};

  Kokkos::parallel_for("fix rigid/nh/small final body",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody) {
      Body &b = d_body(ibody);

      // update vcm by 1/2 step

      double dtfm = l_dtf / b.mass;
      if (scaleflag) {
        b.vcm[0] *= scale_t[0];
        b.vcm[1] *= scale_t[1];
        b.vcm[2] *= scale_t[2];
      }

      b.vcm[0] += dtfm * b.fcm[0];
      b.vcm[1] += dtfm * b.fcm[1];
      b.vcm[2] += dtfm * b.fcm[2];

      // update conjqm, transform to angmom, set velocity again
      // rotational quantities in KK_FLOAT; quat/conjqm/fquat stay double

      KK_FLOAT ex[3] = {(KK_FLOAT)b.ex_space[0],(KK_FLOAT)b.ex_space[1],(KK_FLOAT)b.ex_space[2]};
      KK_FLOAT ey[3] = {(KK_FLOAT)b.ey_space[0],(KK_FLOAT)b.ey_space[1],(KK_FLOAT)b.ey_space[2]};
      KK_FLOAT ez[3] = {(KK_FLOAT)b.ez_space[0],(KK_FLOAT)b.ez_space[1],(KK_FLOAT)b.ez_space[2]};
      KK_FLOAT torque[3] = {(KK_FLOAT)b.torque[0],(KK_FLOAT)b.torque[1],(KK_FLOAT)b.torque[2]};
      KK_FLOAT inertia[3] = {(KK_FLOAT)b.inertia[0],(KK_FLOAT)b.inertia[1],(KK_FLOAT)b.inertia[2]};
      KK_FLOAT tbody[3], mbody[3], angmom[3], omega[3];
      double fquat[4];
      MathExtraKokkos::transpose_matvec(ex,ey,ez,torque,tbody);
      MathExtraKokkos::quatvec(b.quat,tbody,fquat);

      if (scaleflag) {
        b.conjqm[0] = scale_r_l * b.conjqm[0] + dtf2 * fquat[0];
        b.conjqm[1] = scale_r_l * b.conjqm[1] + dtf2 * fquat[1];
        b.conjqm[2] = scale_r_l * b.conjqm[2] + dtf2 * fquat[2];
        b.conjqm[3] = scale_r_l * b.conjqm[3] + dtf2 * fquat[3];
      } else {
        b.conjqm[0] += dtf2 * fquat[0];
        b.conjqm[1] += dtf2 * fquat[1];
        b.conjqm[2] += dtf2 * fquat[2];
        b.conjqm[3] += dtf2 * fquat[3];
      }

      MathExtraKokkos::invquatvec(b.quat,b.conjqm,mbody);
      MathExtraKokkos::matvec(ex,ey,ez,mbody,angmom);

      angmom[0] *= 0.5;
      angmom[1] *= 0.5;
      angmom[2] *= 0.5;

      MathExtraKokkos::angmom_to_omega(angmom,ex,ey,ez,inertia,omega);

      // store the updated rotational state back to the (double) body
      b.angmom[0]=angmom[0]; b.angmom[1]=angmom[1]; b.angmom[2]=angmom[2];
      b.omega[0]=omega[0]; b.omega[1]=omega[1]; b.omega[2]=omega[2];
    });
}

/* ----------------------------------------------------------------------
   accumulate translational and rotational kinetic energies over local bodies
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::nh_akin(double &akin_t_out, double &akin_r_out)
{
  auto d_body = this->d_body;
  double l_akin_t = 0.0, l_akin_r = 0.0;

  Kokkos::parallel_reduce("fix rigid/nh/small akin",
    Range1D(0, nlocal_body),
    KOKKOS_LAMBDA(const int ibody, double &at, double &ar) {
      Body &b = d_body(ibody);
      at += b.mass*(b.vcm[0]*b.vcm[0] + b.vcm[1]*b.vcm[1] + b.vcm[2]*b.vcm[2]);
      ar += b.angmom[0]*b.omega[0] + b.angmom[1]*b.omega[1] + b.angmom[2]*b.omega[2];
    }, l_akin_t, l_akin_r);

  double ke[2],keall[2];
  ke[0] = l_akin_t;
  ke[1] = l_akin_r;
  MPI_Allreduce(ke,keall,2,MPI_DOUBLE,MPI_SUM,world);
  akin_t_out = keall[0];
  akin_r_out = keall[1];
}

/* ----------------------------------------------------------------------
   Nose-Hoover thermostat chain integration (host, global scalar state)
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::nhc_temp_integrate()
{
  if (g_f == 0) return;

  int i,j,k;
  double kt,gfkt_t,gfkt_r,tmp,ms,s,s2;

  kt = boltz * t_target;
  gfkt_t = nf_t * kt;
  gfkt_r = nf_r * kt;

  double t_mass = boltz * t_target / (t_freq * t_freq);
  q_t[0] = nf_t * t_mass;
  q_r[0] = nf_r * t_mass;
  for (i = 1; i < t_chain; i++)
    q_t[i] = q_r[i] = t_mass;

  f_eta_t[0] = (akin_t * mvv2e - gfkt_t) / q_t[0];
  f_eta_r[0] = (akin_r * mvv2e - gfkt_r) / q_r[0];

  for (i = 0; i < t_iter; i++) {
    for (j = 0; j < t_order; j++) {

      eta_dot_t[t_chain-1] += wdti2[j] * f_eta_t[t_chain-1];
      eta_dot_r[t_chain-1] += wdti2[j] * f_eta_r[t_chain-1];

      for (k = 1; k < t_chain; k++) {
        tmp = wdti4[j] * eta_dot_t[t_chain-k];
        ms = maclaurin_series(tmp);
        s = exp(-1.0 * tmp);
        s2 = s * s;
        eta_dot_t[t_chain-k-1] = eta_dot_t[t_chain-k-1] * s2 +
          wdti2[j] * f_eta_t[t_chain-k-1] * s * ms;

        tmp = wdti4[j] * eta_dot_r[t_chain-k];
        ms = maclaurin_series(tmp);
        s = exp(-1.0 * tmp);
        s2 = s * s;
        eta_dot_r[t_chain-k-1] = eta_dot_r[t_chain-k-1] * s2 +
          wdti2[j] * f_eta_r[t_chain-k-1] * s * ms;
      }

      for (k = 0; k < t_chain; k++) {
        eta_t[k] += wdti1[j] * eta_dot_t[k];
        eta_r[k] += wdti1[j] * eta_dot_r[k];
      }

      for (k = 1; k < t_chain; k++) {
        f_eta_t[k] = q_t[k-1] * eta_dot_t[k-1] * eta_dot_t[k-1] - kt;
        f_eta_t[k] /= q_t[k];
        f_eta_r[k] = q_r[k-1] * eta_dot_r[k-1] * eta_dot_r[k-1] - kt;
        f_eta_r[k] /= q_r[k];
      }

      for (k = 0; k < t_chain-1; k++) {
        tmp = wdti4[j] * eta_dot_t[k+1];
        ms = maclaurin_series(tmp);
        s = exp(-1.0 * tmp);
        s2 = s * s;
        eta_dot_t[k] = eta_dot_t[k] * s2 + wdti2[j] * f_eta_t[k] * s * ms;
        tmp = q_t[k] * eta_dot_t[k] * eta_dot_t[k] - kt;
        f_eta_t[k+1] = tmp / q_t[k+1];

        tmp = wdti4[j] * eta_dot_r[k+1];
        ms = maclaurin_series(tmp);
        s = exp(-1.0 * tmp);
        s2 = s * s;
        eta_dot_r[k] = eta_dot_r[k] * s2 + wdti2[j] * f_eta_r[k] * s * ms;
        tmp = q_r[k] * eta_dot_r[k] * eta_dot_r[k] - kt;
        f_eta_r[k+1] = tmp / q_r[k+1];
      }

      eta_dot_t[t_chain-1] += wdti2[j] * f_eta_t[t_chain-1];
      eta_dot_r[t_chain-1] += wdti2[j] * f_eta_r[t_chain-1];
    }
  }
}

/* ----------------------------------------------------------------------
   Nose-Hoover barostat chain integration (host, global scalar state)
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::nhc_press_integrate()
{
  int i,j,k;
  double tmp,s,s2,ms,kecurrent;
  double kt = boltz * t_target;
  double lkt_press = kt;

  int dimension = domain->dimension;
  double tb_mass = kt / (p_freq_max * p_freq_max);
  q_b[0] = dimension * dimension * tb_mass;
  for (i = 1; i < p_chain; i++) {
    q_b[i] = tb_mass;
    f_eta_b[i] = q_b[i-1] * eta_dot_b[i-1] * eta_dot_b[i-1] - kt;
    f_eta_b[i] /= q_b[i];
  }

  kecurrent = 0.0;
  for (i = 0; i < 3; i++)
    if (p_flag[i]) {
      epsilon_mass[i] = (g_f + dimension) * kt / (p_freq[i] * p_freq[i]);
      kecurrent += epsilon_mass[i] * epsilon_dot[i] * epsilon_dot[i];
    }
  kecurrent /= pdim;

  f_eta_b[0] = (kecurrent - lkt_press) / q_b[0];

  for (i = 0; i < t_iter; i++) {
    for (j = 0; j < t_order; j++) {

      eta_dot_b[p_chain-1] += wdti2[j] * f_eta_b[p_chain-1];

      for (k = 1; k < p_chain; k++) {
        tmp = wdti4[j] * eta_dot_b[p_chain-k];
        ms = maclaurin_series(tmp);
        s = exp(-0.5 * tmp);
        s2 = s * s;
        eta_dot_b[p_chain-k-1] = eta_dot_b[p_chain-k-1] * s2 +
          wdti2[j] * f_eta_b[p_chain-k-1] * s * ms;
      }

      for (k = 0; k < p_chain; k++)
        eta_b[k] += wdti1[j] * eta_dot_b[k];

      for (k = 1; k < p_chain; k++) {
        f_eta_b[k] = q_b[k-1] * eta_dot_b[k-1] * eta_dot_b[k-1] - kt;
        f_eta_b[k] /= q_b[k];
      }

      for (k = 0; k < p_chain-1; k++) {
        tmp = wdti4[j] * eta_dot_b[k+1];
        ms = maclaurin_series(tmp);
        s = exp(-0.5 * tmp);
        s2 = s * s;
        eta_dot_b[k] = eta_dot_b[k] * s2 + wdti2[j] * f_eta_b[k] * s * ms;
        tmp = q_b[k] * eta_dot_b[k] * eta_dot_b[k] - kt;
        f_eta_b[k+1] = tmp / q_b[k+1];
      }

      eta_dot_b[p_chain-1] += wdti2[j] * f_eta_b[p_chain-1];
    }
  }
}

/* ----------------------------------------------------------------------
   compute kinetic energy in the extended Hamiltonian
   conserved quantity = sum of returned energy and potential energy
-----------------------------------------------------------------------*/

template<class DeviceType>
double FixRigidNHSmallKokkos<DeviceType>::compute_scalar()
{
  int i,k;
  double kt = boltz * t_target;
  double energy,ke_t,ke_q,tmp,Pkq[4];

  double *vcm,*quat;

  // body data lives on the device during a run; flush to host
  copy_body_host();

  ke_t = 0.0;
  ke_q = 0.0;

  for (i = 0; i < nlocal_body; i++) {
    vcm = body[i].vcm;
    quat = body[i].quat;
    ke_t += body[i].mass * (vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);

    for (k = 1; k < 4; k++) {
      if (k == 1) {
        Pkq[0] = -quat[1];
        Pkq[1] =  quat[0];
        Pkq[2] =  quat[3];
        Pkq[3] = -quat[2];
      } else if (k == 2) {
        Pkq[0] = -quat[2];
        Pkq[1] = -quat[3];
        Pkq[2] =  quat[0];
        Pkq[3] =  quat[1];
      } else if (k == 3) {
        Pkq[0] = -quat[3];
        Pkq[1] =  quat[2];
        Pkq[2] = -quat[1];
        Pkq[3] =  quat[0];
      }

      tmp = body[i].conjqm[0]*Pkq[0] + body[i].conjqm[1]*Pkq[1] +
        body[i].conjqm[2]*Pkq[2] + body[i].conjqm[3]*Pkq[3];
      tmp *= tmp;

      if (fabs(body[i].inertia[k-1]) < 1e-6) tmp = 0.0;
      else tmp /= (8.0 * body[i].inertia[k-1]);
      ke_q += tmp;
    }
  }

  double ke[2],keall[2];
  ke[0] = ke_t;
  ke[1] = ke_q;
  MPI_Allreduce(ke,keall,2,MPI_DOUBLE,MPI_SUM,world);
  ke_t = keall[0];
  ke_q = keall[1];

  energy = (ke_t + ke_q) * mvv2e;

  if (tstat_flag) {
    energy += kt * (nf_t * eta_t[0] + nf_r * eta_r[0]);

    for (i = 1; i < t_chain; i++)
      energy += kt * (eta_t[i] + eta_r[i]);

    for (i = 0;  i < t_chain; i++) {
      energy += 0.5 * q_t[i] * (eta_dot_t[i] * eta_dot_t[i]);
      energy += 0.5 * q_r[i] * (eta_dot_r[i] * eta_dot_r[i]);
    }
  }

  if (pstat_flag) {
    double e = 0.0;
    for (i = 0; i < 3; i++)
      if (p_flag[i])
        e += epsilon_mass[i] * epsilon_dot[i] * epsilon_dot[i];
    energy += e*(0.5/pdim);

    double vol;
    if (domain->dimension == 2) vol = domain->xprd * domain->yprd;
    else vol = domain->xprd * domain->yprd * domain->zprd;

    double p0 = (p_target[0] + p_target[1] + p_target[2]) / 3.0;
    energy += p0 * vol / nktv2p;

    for (i = 0;  i < p_chain; i++) {
      energy += kt * eta_b[i];
      energy += 0.5 * q_b[i] * (eta_dot_b[i] * eta_dot_b[i]);
    }
  }

  return energy;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::couple()
{
  double *tensor = pressure->vector;

  if (pstyle == ISO) {
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
  } else if (pcouple == XYZ) {
    double ave = 1.0/3.0 * (tensor[0] + tensor[1] + tensor[2]);
    p_current[0] = p_current[1] = p_current[2] = ave;
  } else if (pcouple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
  } else if (pcouple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
  } else if (pcouple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
  } else {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::remap()
{
  int i;
  double oldlo,oldhi,ctr,expfac;

  // remap dilates atom positions on the host; flush x to host and mark dirty
  atomKK->sync(Host, X_MASK);

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < 3; i++) epsilon[i] += dtq * epsilon_dot[i];

  if (allremap) domain->x2lamda(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->x2lamda(x[i],x[i]);
  }

  for (auto &ifix : rfix) ifix->deform(0);

  for (i = 0; i < 3; i++) {
    if (p_flag[i]) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      expfac = exp(dtq * epsilon_dot[i]);
      domain->boxlo[i] = (oldlo-ctr)*expfac + ctr;
      domain->boxhi[i] = (oldhi-ctr)*expfac + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();

  if (allremap) domain->lamda2x(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->lamda2x(x[i],x[i]);
  }

  for (auto &ifix : rfix) ifix->deform(1);

  atomKK->modified(Host, X_MASK);
}

/* ----------------------------------------------------------------------
   compute target temperature and kinetic energy
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::compute_temp_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  t_target = t_start + delta * (t_stop-t_start);
}

/* ----------------------------------------------------------------------
   compute hydrostatic target pressure
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::compute_press_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  p_hydro = 0.0;
  for (int i = 0; i < 3; i++)
    if (p_flag[i]) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      p_hydro += p_target[i];
    }
  p_hydro /= pdim;
}

/* ----------------------------------------------------------------------
   update epsilon_dot
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::nh_epsilon_dot()
{
  if (g_f == 0) return;

  int i;
  double volume,scale,f_epsilon;

  if (domain->dimension == 2) volume = domain->xprd*domain->yprd;
  else volume = domain->xprd*domain->yprd*domain->zprd;

  mtk_term1 = (akin_t + akin_r) * mvv2e / g_f;

  scale = exp(-1.0 * dtq * eta_dot_b[0]);

  for (i = 0; i < 3; i++)
    if (p_flag[i]) {
      f_epsilon = (p_current[i]-p_hydro)*volume / nktv2p + mtk_term1;
      f_epsilon /= epsilon_mass[i];
      epsilon_dot[i] += dtq * f_epsilon;
      epsilon_dot[i] *= scale;
    }

  mtk_term2 = 0.0;
  for (i = 0; i < 3; i++)
    if (p_flag[i]) mtk_term2 += epsilon_dot[i];
  mtk_term2 /= g_f;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::compute_dof()
{
  // reads host body[].inertia; caller ensures body[] is valid on host
  int dimension = domain->dimension;

  nf_t = dimension * nlocal_body;
  if (dimension == 3) {
    nf_r = dimension * nlocal_body;
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      Body *b = &body[ibody];
      for (int k = 0; k < dimension; k++)
        if (fabs(b->inertia[k]) < EPSILON) nf_r--;
    }
  } else if (dimension == 2) nf_r = nlocal_body;

  double nf[2], nfall[2];
  nf[0] = nf_t;
  nf[1] = nf_r;
  MPI_Allreduce(nf,nfall,2,MPI_DOUBLE,MPI_SUM,world);
  nf_t = (int)nfall[0];
  nf_r = (int)nfall[1];

  g_f = nf_t + nf_r;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::write_restart(FILE *fp)
{
  if (tstat_flag == 0 && pstat_flag == 0) return;

  int nsize = 2; // tstat_flag and pstat_flag

  if (tstat_flag) {
    nsize += 1;         // t_chain
    nsize += 4*t_chain; // eta_t, eta_r, eta_dot_t, eta_dot_r
  }

  if (pstat_flag) {
    nsize += 7;         // p_chain, epsilon(3) and epsilon_dot(3)
    nsize += 2*p_chain;
  }

  double *list;
  memory->create(list,nsize,"rigid_nh:list");

  int n = 0;

  list[n++] = tstat_flag;
  if (tstat_flag) {
    list[n++] = t_chain;
    for (int i = 0; i < t_chain; i++) {
      list[n++] = eta_t[i];
      list[n++] = eta_r[i];
      list[n++] = eta_dot_t[i];
      list[n++] = eta_dot_r[i];
    }
  }

  list[n++] = pstat_flag;
  if (pstat_flag) {
    list[n++] = epsilon[0];
    list[n++] = epsilon[1];
    list[n++] = epsilon[2];
    list[n++] = epsilon_dot[0];
    list[n++] = epsilon_dot[1];
    list[n++] = epsilon_dot[2];

    list[n++] = p_chain;
    for (int i = 0; i < p_chain; i++) {
      list[n++] = eta_b[i];
      list[n++] = eta_dot_b[i];
    }
  }

  if (comm->me == 0) {
    int size = nsize*sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),nsize,fp);
  }

  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::restart(char *buf)
{
  int n = 0;
  auto *list = (double *) buf;
  int flag = static_cast<int> (list[n++]);

  if (flag) {
    int m = static_cast<int> (list[n++]);
    if (tstat_flag && m == t_chain) {
      for (int i = 0; i < t_chain; i++) {
        eta_t[i] = list[n++];
        eta_r[i] = list[n++];
        eta_dot_t[i] = list[n++];
        eta_dot_r[i] = list[n++];
      }
    } else n += 4*m;
  }

  flag = static_cast<int> (list[n++]);
  if (flag) {
    epsilon[0] = list[n++];
    epsilon[1] = list[n++];
    epsilon[2] = list[n++];
    epsilon_dot[0] = list[n++];
    epsilon_dot[1] = list[n++];
    epsilon_dot[2] = list[n++];

    int m = static_cast<int> (list[n++]);
    if (pstat_flag && m == p_chain) {
      for (int i = 0; i < p_chain; i++) {
        eta_b[i] = list[n++];
        eta_dot_b[i] = list[n++];
      }
    } else n += 2*m;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidNHSmallKokkos<DeviceType>::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tcomputeflag) {
      modify->delete_compute(id_temp);
      tcomputeflag = 0;
    }
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);

    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature) error->all(FLERR,"Could not find fix_modify temperature ID {}", id_temp);

    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify temperature ID {} does not compute temperature", id_temp);
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for fix modify is not for group all");

    if (pstat_flag) {
      pressure = modify->get_compute_by_id(id_press);
      if (!pressure) error->all(FLERR,"Pressure ID {} for fix modify does not exist", id_press);
      pressure->reset_extra_compute_fix(id_temp);
    }

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (!pstat_flag) error->all(FLERR,"Illegal fix_modify command");
    if (pcomputeflag) {
      modify->delete_compute(id_press);
      pcomputeflag = 0;
    }
    delete[] id_press;
    id_press = utils::strdup(arg[1]);

    pressure = modify->get_compute_by_id(id_press);
    if (!pressure) error->all(FLERR,"Could not find fix_modify pressure ID {}", id_press);

    if (pressure->pressflag == 0)
      error->all(FLERR,"Fix_modify pressure ID {} does not compute pressure", id_press);
    return 2;
  }

  return FixRigidSmall::modify_param(narg,arg);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::allocate_chain()
{
  if (tstat_flag) {
    q_t = new double[t_chain];
    q_r = new double[t_chain];
    eta_t = new double[t_chain];
    eta_r = new double[t_chain];
    eta_dot_t = new double[t_chain];
    eta_dot_r = new double[t_chain];
    f_eta_t = new double[t_chain];
    f_eta_r = new double[t_chain];
  }

  if (pstat_flag) {
    q_b = new double[p_chain];
    eta_b = new double[p_chain];
    eta_dot_b = new double[p_chain];
    f_eta_b = new double[p_chain];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::allocate_order()
{
  w = new double[t_order];
  wdti1 = new double[t_order];
  wdti2 = new double[t_order];
  wdti4 = new double[t_order];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::deallocate_chain()
{
  if (tstat_flag) {
    delete[] q_t;
    delete[] q_r;
    delete[] eta_t;
    delete[] eta_r;
    delete[] eta_dot_t;
    delete[] eta_dot_r;
    delete[] f_eta_t;
    delete[] f_eta_r;
  }

  if (pstat_flag) {
    delete[] q_b;
    delete[] eta_b;
    delete[] eta_dot_b;
    delete[] f_eta_b;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidNHSmallKokkos<DeviceType>::deallocate_order()
{
  delete[] w;
  delete[] wdti1;
  delete[] wdti2;
  delete[] wdti4;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidNHSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidNHSmallKokkos<LMPHostType>;
#endif
}
