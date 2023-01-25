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
FixStyle(mdi/qmmm,FixMDIQMMM);
// clang-format on
#else

#ifndef LMP_FIX_MDI_QMMM_H
#define LMP_FIX_MDI_QMMM_H

#include "fix.h"
#include <mdi.h>

namespace LAMMPS_NS {

class FixMDIQMMM : public Fix {
 public:
  FixMDIQMMM(class LAMMPS *, int, char **);
  virtual ~FixMDIQMMM();
  int setmask() override;

  void init() override;
  void setup(int) override;
  void setup_post_neighbor() override;
  void setup_pre_force(int) override;
  void post_neighbor() override;
  void pre_force(int) override;
  void post_force(int) override;
  void min_setup(int) override;
  void min_post_neighbor() override;
  void min_pre_force(int) override;
  void min_post_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  double compute_scalar() override;
  double compute_vector(int) override;
  double memory_usage() override;
  
 private:
  int nprocs;
  int virialflag, connectflag;
  int plugin;
  int maxlocal;
  int sumflag;
  int *elements;
  int mode;            // QMMM method = DIRECT or POTENTIAL

  int lmpunits;        // REAL, METAL, or NATIVE
  int first_send;      // 1 until initial info passed to MDI engine
  
  double qm_energy;
  double qm_virial[9], qm_virial_symmetric[6];

  MDI_Comm mdicomm;

  class Pair *pair_coul;    // ptr to instance of pair coul variant
  
  // data for QM portion

  int nqm;                   // # of QM atoms
  tagint *qmIDs;             // IDs of QM atoms in ascending order
  double **xqm,**fqm;        // QM coords and forces
  double *qqm;               // QM charges
  int *eqm;                  // QM atom atomic numbers
  double *qpotential;        // Coulomb potential
  double **xqm_mine;         // same values for QM atoms I own
  double *qqm_mine;
  double *qpotential_mine;
  int *eqm_mine;
  int *qm2owned;             // index of local atom for each QM atom
                             // index = -1 if this proc does not own
  
  double *ecoul;             // peratom Coulombic energy from LAMMPS
  int ncoulmax;              // length of ecoul

  // data for MM portion

  int nmm;                   // # of MM atoms
  tagint *mmIDs;             // IDs of MM atoms in ascending order
  double **xmm,**fmm;        // MM coords and forces
  double *qmm;               // MM charges
  int *emm;                  // MM atom atomic numbers
  double **xmm_mine;         // same values for MM atoms I own
  double *qmm_mine;
  int *emm_mine;
  int *mm2owned;             // index of local atom for each MM atom
                             // index = -1 if this proc does not own

  // unit conversion factors

  double lmp2mdi_length, mdi2lmp_length;
  double lmp2mdi_energy, mdi2lmp_energy;
  double lmp2mdi_force, mdi2lmp_force;
  double lmp2mdi_pressure, mdi2lmp_pressure;

  // methods

  void post_force_direct(int);
  void post_force_potential(int);

  void create_qm_list();
  void create_mm_list();
  
  void set_qm2owned();
  void set_mm2owned();
  
  void set_eqm();
  void set_tqm();
  void set_qqm();
  void set_xqm();

  void set_emm();
  void set_qmm();
  void set_xmm();
  
  void send_box();
  void unit_conversions();

};

}    // namespace LAMMPS_NS

#endif
#endif
