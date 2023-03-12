/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  virtual ~FixMDIQM();
  int setmask() override;

  void init() override;
  void setup(int) override;
  void setup_post_neighbor() override;

  void post_neighbor() override;
  void post_force(int) override;

  void min_setup(int) override;
  void min_post_neighbor() override;
  void min_post_force(int) override;

  double compute_scalar() override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  int nprocs;
  int every, virialflag, addflag, connectflag;
  int plugin;
  int sumflag,mcflag;
  char *id_mcfix;
  int *mc_active_ptr,*exclusion_group_ptr;
  int *elements;

  double qm_cell[9],qm_cell_displ[3];

  double qm_energy;
  double qm_virial[9], qm_virial_symmetric[6];

  MDI_Comm mdicomm;
  int natoms_exists,celldispl_exists,elements_exists,types_exists;
  int stress_exists;

  int nmax;

  // unit conversion factors

  int lmpunits;

  double lmp2mdi_length, mdi2lmp_length;
  double lmp2mdi_energy, mdi2lmp_energy;
  double lmp2mdi_force, mdi2lmp_force;
  double lmp2mdi_pressure, mdi2lmp_pressure;

  // QM atom data structs

  int nqm,nqm_last,max_nqm;
  int nexclude;

  tagint *qmIDs;
  int *qm2owned;

  int *eqm,*eqm_mine;
  int *tqm,*tqm_mine;
  double **xqm,**xqm_mine;
  double **fqm;

  // methods

  void reallocate();

  int set_nqm();
  void create_qm_list();
  void set_qm2owned();
  void set_box();
  void set_xqm();
  void set_eqm();
  void set_tqm();

  void send_natoms();
  void send_types();
  void send_elements();
  void send_box();

  void unit_conversions();
};

}    // namespace LAMMPS_NS

#endif
#endif
