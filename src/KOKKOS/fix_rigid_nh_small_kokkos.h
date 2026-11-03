/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_RIGID_NH_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_NH_SMALL_KOKKOS_H

#include "fix_rigid_small_kokkos.h"

namespace LAMMPS_NS {

// Kokkos port of the Nose-Hoover (thermostat/barostat) rigid-body integrator.
//
// The per-body velocity-Verlet update (with the no_squish quaternion scheme and
// the thermostat/barostat scale factors) runs as device kernels over the body
// DualView owned by FixRigidSmallKokkos.  The Nose-Hoover chain integration and
// pressure/temperature coupling are global scalar operations and stay on the
// host (they touch only fixed-size 1d chain state, not per-atom/per-body data),
// matching the design used by fix_nh_kokkos.

template <class DeviceType> class FixRigidNHSmallKokkos : public FixRigidSmallKokkos<DeviceType> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixRigidNHSmallKokkos(class LAMMPS *, int, char **);
  ~FixRigidNHSmallKokkos() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  double compute_scalar() override;
  int modify_param(int, char **) override;
  void write_restart(FILE *) override;
  void restart(char *buf) override;
  void reset_target(double) override;

  // run the per-body initial/final Nose-Hoover update on the device
  // (must be public: nvcc's extended __host__ __device__ lambda extension
  // requires the enclosing function to have public access)
  void nh_initial_integrate_bodies(const double scale_r, const double scale_t[3],
                                   const double scale_v[3]);
  void nh_final_integrate_bodies(const double scale_r, const double scale_t[3]);
  void nh_akin(double &akin_t_out, double &akin_r_out);

 protected:
  using Range1D = Kokkos::RangePolicy<DeviceType>;

  // bring inherited (dependent base) members into scope so the host scalar
  // chain machinery copied from FixRigidNHSmall compiles unqualified
  using FixRigidSmallKokkos<DeviceType>::atom;
  using FixRigidSmallKokkos<DeviceType>::atomKK;
  using FixRigidSmallKokkos<DeviceType>::domain;
  using FixRigidSmallKokkos<DeviceType>::force;
  using FixRigidSmallKokkos<DeviceType>::update;
  using FixRigidSmallKokkos<DeviceType>::modify;
  using FixRigidSmallKokkos<DeviceType>::group;
  using FixRigidSmallKokkos<DeviceType>::comm;
  using FixRigidSmallKokkos<DeviceType>::error;
  using FixRigidSmallKokkos<DeviceType>::memory;
  using FixRigidSmallKokkos<DeviceType>::world;
  using FixRigidSmallKokkos<DeviceType>::execution_space;
  using FixRigidSmallKokkos<DeviceType>::datamask_read;
  using FixRigidSmallKokkos<DeviceType>::datamask_modify;

  using FixRigidSmallKokkos<DeviceType>::style;
  using FixRigidSmallKokkos<DeviceType>::dynamic;
  using FixRigidSmallKokkos<DeviceType>::virial;
  using FixRigidSmallKokkos<DeviceType>::evflag;

  using FixRigidSmallKokkos<DeviceType>::body;
  using FixRigidSmallKokkos<DeviceType>::nlocal_body;
  using FixRigidSmallKokkos<DeviceType>::nghost_body;
  using FixRigidSmallKokkos<DeviceType>::nmax_body;
  using FixRigidSmallKokkos<DeviceType>::bodysize;
  using FixRigidSmallKokkos<DeviceType>::commflag;
  using FixRigidSmallKokkos<DeviceType>::tstat_flag;
  using FixRigidSmallKokkos<DeviceType>::pstat_flag;
  using FixRigidSmallKokkos<DeviceType>::t_chain;
  using FixRigidSmallKokkos<DeviceType>::t_iter;
  using FixRigidSmallKokkos<DeviceType>::t_order;
  using FixRigidSmallKokkos<DeviceType>::p_chain;
  using FixRigidSmallKokkos<DeviceType>::p_flag;
  using FixRigidSmallKokkos<DeviceType>::p_freq;
  using FixRigidSmallKokkos<DeviceType>::p_start;
  using FixRigidSmallKokkos<DeviceType>::p_stop;
  using FixRigidSmallKokkos<DeviceType>::p_period;
  using FixRigidSmallKokkos<DeviceType>::p_current;
  using FixRigidSmallKokkos<DeviceType>::p_target;
  using FixRigidSmallKokkos<DeviceType>::pstyle;
  using FixRigidSmallKokkos<DeviceType>::pcouple;
  using FixRigidSmallKokkos<DeviceType>::t_start;
  using FixRigidSmallKokkos<DeviceType>::t_stop;
  using FixRigidSmallKokkos<DeviceType>::t_period;
  using FixRigidSmallKokkos<DeviceType>::earlyflag;
  using FixRigidSmallKokkos<DeviceType>::langflag;
  using FixRigidSmallKokkos<DeviceType>::dtv;
  using FixRigidSmallKokkos<DeviceType>::dtf;
  using FixRigidSmallKokkos<DeviceType>::dtq;
  using FixRigidSmallKokkos<DeviceType>::allremap;
  using FixRigidSmallKokkos<DeviceType>::dilate_group_bit;
  using FixRigidSmallKokkos<DeviceType>::id_dilate;
  using FixRigidSmallKokkos<DeviceType>::setupflag;

  // device machinery from FixRigidSmallKokkos
  using FixRigidSmallKokkos<DeviceType>::commKK;
  using FixRigidSmallKokkos<DeviceType>::d_body;
  using FixRigidSmallKokkos<DeviceType>::k_body;
  using FixRigidSmallKokkos<DeviceType>::set_xv_kokkos;
  using FixRigidSmallKokkos<DeviceType>::setup_device_push;
  using FixRigidSmallKokkos<DeviceType>::compute_forces_and_torques_kokkos;
  using FixRigidSmallKokkos<DeviceType>::copy_body_host;
  using FixRigidSmallKokkos<DeviceType>::copy_body_device;

  typedef typename FixRigidSmall::Body Body;

  double boltz, nktv2p, mvv2e;    // boltzman constant, conversion factors

  int nf_t, nf_r;                       // trans/rot degrees of freedom
  double *w, *wdti1, *wdti2, *wdti4;    // Yoshida-Suzuki coefficients
  double *q_t, *q_r;                    // trans/rot thermostat masses
  double *eta_t, *eta_r;                // trans/rot thermostat positions
  double *eta_dot_t, *eta_dot_r;        // trans/rot thermostat velocities
  double *f_eta_t, *f_eta_r;            // trans/rot thermostat forces

  double epsilon_mass[3], *q_b;         // baro/thermo masses
  double epsilon[3], *eta_b;            // baro/thermo positions
  double epsilon_dot[3], *eta_dot_b;    // baro/thermo velocities
  double *f_eta_b;                      // thermo forces
  double akin_t, akin_r;                // translational/rotational kinetic energies

  int kspace_flag;            // 1 if KSpace invoked, 0 if not
  std::vector<Fix *> rfix;    // indices of rigid fixes

  double vol0;          // reference volume
  double t0;            // reference temperature
  int pdim, g_f;        // number of barostatted dims, total DoFs
  double p_hydro;       // hydrostatic target pressure
  double p_freq_max;    // maximum barostat frequency

  double mtk_term1, mtk_term2;    // Martyna-Tobias-Klein corrections

  double t_target, t_current;
  double t_freq;

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tcomputeflag, pcomputeflag;    // 1 = compute was created by fix. 0 = external

  void couple();
  void remap();
  void nhc_temp_integrate();
  void nhc_press_integrate();

  virtual void compute_temp_target();
  void compute_press_target();
  void nh_epsilon_dot();
  void compute_dof();

  void allocate_chain();
  void allocate_order();
  void deallocate_chain();
  void deallocate_order();

  inline double maclaurin_series(double);
};

template <class DeviceType>
inline double FixRigidNHSmallKokkos<DeviceType>::maclaurin_series(double x)
{
  double x2, x4;
  x2 = x * x;
  x4 = x2 * x2;
  return (1.0 + (1.0 / 6.0) * x2 + (1.0 / 120.0) * x4 + (1.0 / 5040.0) * x2 * x4 +
          (1.0 / 362880.0) * x4 * x4);
}

}    // namespace LAMMPS_NS

#endif
