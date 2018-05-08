/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------
   USER-BOCS written by: Nicholas J. H. Dunn and Michael R. DeLyser
   from The Pennsylvania State University
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(PRESSURE/BOCS,ComputePressureBocs)

#else


#ifndef LMP_COMPUTE_PRESSURE_BOCS_H
#define LMP_COMPUTE_PRESSURE_BOCS_H

#include "compute.h"

namespace LAMMPS_NS {
// ComputePressure -> ComputePressureBocs MRD NJD
class ComputePressureBocs : public Compute {
 public:
  ComputePressureBocs(class LAMMPS *, int, char **);
  virtual ~ComputePressureBocs();
  virtual void init();
  virtual double compute_scalar();
  virtual void compute_vector();
  void reset_extra_compute_fix(const char *);

  double compute_cg_scalar();
  double get_cg_p_corr(int, double *, int, double, double);
  double get_cg_fluct(double, double);
  void send_cg_info(int, int, double*, int, double);
  void send_cg_info(int, double **, int);
  double get_cg_p_corr(double **, int, double);
  double find_index(double*  , double);

 protected:
  double boltz,nktv2p,inv_volume;
  int nvirial,dimension;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int fixflag,kspaceflag;

// NJD MRD
  int p_basis_type;
  int p_match_flag;
  double vavg;
  int N_mol;
  int N_basis;
  double *phi_coeff;
  double ** splines;
  int spline_length;

  void virial_compute(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pressure must use group all

Virial contributions computed by potentials (pair, bond, etc) are
computed on all atoms.

E: Could not find compute pressure temperature ID

The compute ID for calculating temperature does not exist.

E: Compute pressure temperature ID does not compute temperature

The compute ID assigned to a pressure computation must compute
temperature.

E: Compute pressure requires temperature ID to include kinetic energy

The keflag cannot be used unless a temperature compute is provided.

E: Virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

E: Must use 'kspace_modify pressure/scalar no' for tensor components with kspace_style msm

Otherwise MSM will compute only a scalar pressure.  See the kspace_modify
command for details on this setting.

*/
