/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(mdi/qm,FixMDIQM);
// clang-format on
#else

#ifndef LMP_FIX_MDI_QM_H
#define LMP_FIX_MDI_QM_H

#include "fix.h"
#include <mdi.h>

namespace LAMMPS_NS {

class FixMDIQM : public Fix {
 public:
  FixMDIQM(class LAMMPS *, int, char **);
  ~FixMDIQM();
  int setmask();

  void init();
  void setup(int);
  void post_force(int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  int nprocs;
  int every, virialflag, addflag, connectflag;
  int plugin;
  int maxlocal;
  int sumflag;
  int *elements;

  double qm_energy;
  int lmpunits;
  double qm_virial[9], qm_virial_symmetric[6];
  double **fqm;

  MDI_Comm mdicomm;

  // unit conversion factors

  double lmp2mdi_length, mdi2lmp_length;
  double lmp2mdi_energy, mdi2lmp_energy;
  double lmp2mdi_force, mdi2lmp_force;
  double lmp2mdi_pressure, mdi2lmp_pressure;
  double lmp2mdi_velocity, mdi2lmp_velocity;

  // buffers for MDI comm

  int maxbuf;
  int *ibuf1, *ibuf1all;
  double *buf3, *buf3all;

  // methods

  void reallocate();
  void send_types();
  void send_elements();
  void send_box();
  void unit_conversions();
};

}    // namespace LAMMPS_NS

#endif
#endif
