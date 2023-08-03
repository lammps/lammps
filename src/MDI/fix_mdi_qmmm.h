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
  int mode;    // QMMM method = DIRECT or POTENTIAL

  double qm_cell[9], qm_cell_displ[3];

  double qm_energy;
  double qm_virial[9], qm_virial_symmetric[6];

  MDI_Comm mdicomm;
  int natoms_exists, celldispl_exists, elements_exists, types_exists;
  int stress_exists, pe_exists, keelec_exists;

  int nmax;

  class Pair *pair_coul;    // ptr to instance of pair coul variant

  // QM atom data structs

  int nqm;    // # of QM atoms
  int nqm_last, max_nqm;

  tagint *qmIDs;    // IDs of QM atoms in ascending order
  int *qm2owned;    // index of local atom for each QM atom
                    // index = -1 if this proc does not own

  int *eqm, *eqm_mine;
  int *tqm, *tqm_mine;
  double **xqm, **xqm_mine;
  double *qqm, *qqm_mine;
  double *qpotential, *qpotential_mine;
  double **fqm;

  double *ecoul;    // peratom Coulombic energy from LAMMPS
  int ncoulmax;     // length of ecoul

  // MM atom data structs

  int nmm;    // # of MM atoms
  int nmm_last, max_nmm;

  tagint *mmIDs;    // IDs of MM atoms in ascending order
  int *mm2owned;    // index of local atom for each MM atom
                    // index = -1 if this proc does not own

  int *emm, *emm_mine;
  int *tmm, *tmm_mine;
  double **xmm, **xmm_mine;
  double *qmm, *qmm_mine;
  double **fmm;

  // unit conversion factors

  int lmpunits;    // REAL, METAL, or NATIVE

  double lmp2mdi_length, mdi2lmp_length;
  double lmp2mdi_energy, mdi2lmp_energy;
  double lmp2mdi_force, mdi2lmp_force;
  double lmp2mdi_pressure, mdi2lmp_pressure;

  // methods

  void post_force_direct(int);
  void post_force_potential(int);

  void reallocate_qm();
  void reallocate_mm();

  int set_nqm();
  int set_nmm();

  void create_qm_list();
  void create_mm_list();

  void set_qm2owned();
  void set_mm2owned();

  void set_box();

  void set_xqm(int);
  void set_eqm();
  void set_tqm();
  void set_qqm();

  void set_xmm();
  void set_emm();
  void set_tmm();
  void set_qmm();

  void send_natoms_qm();
  void send_types_qm();
  void send_elements_qm();

  void send_box();
  void send_natoms_mm();
  void send_types_mm();
  void send_elements_mm();
  void send_charges_mm();

  void request_qm_energy();

  void unit_conversions();
};

}    // namespace LAMMPS_NS

#endif
#endif
