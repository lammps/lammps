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

#ifdef FIX_CLASS
// clang-format off
FixStyle(langevin,FixLangevin);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_H
#define LMP_FIX_LANGEVIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevin : public Fix {
 public:
  FixLangevin(class LAMMPS *, int, char **);
  ~FixLangevin() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void end_of_step() override;
  void reset_target(double) override;
  void reset_dt() override;
  int modify_param(int, char **) override;
  double compute_scalar() override;
  double memory_usage() override;
  void *extract(const char *, int &) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 protected:
  int gjfflag, nvalues, osflag, oflag, tallyflag, zeroflag, tbiasflag;
  int flangevin_allocated;
  double ascale;
  double t_start, t_stop, t_period, t_target;
  double *gfactor1, *gfactor2, *ratio;
  double energy, energy_onestep;
  double tsqrt;
  int tstyle, tvar;
  double gjfa, gjfsib;    //gjf a and gjf sqrt inverse b
  char *tstr;

  class AtomVecEllipsoid *avec;

  int maxatom1, maxatom2;
  double **flangevin;
  double *tforce;
  double **franprev;
  double **lv;    //half step velocity

  char *id_temp;
  class Compute *temperature;

  int nlevels_respa;
  class RanMars *random;
  int seed;

  template <int Tp_TSTYLEATOM, int Tp_GJF, int Tp_TALLY, int Tp_BIAS, int Tp_RMASS, int Tp_ZERO>
  void post_force_templated();

  void omega_thermostat();
  void angmom_thermostat();
  void compute_target();
};

}    // namespace LAMMPS_NS

#endif
#endif
