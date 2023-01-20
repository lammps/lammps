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
  int mode;            // DIRECT or POTENTIAL

  int qflag;           // 1 if per-atom charge defined, 0 if not
  int qm_init;         // 1 if QM code and qqm are initialized, 0 if not

  double qm_energy;
  int lmpunits;
  double qm_virial[9], qm_virial_symmetric[6];

  MDI_Comm mdicomm;

  class Pair *pair_coul;    // ptr to instance of pair coul variant
  
  // data for QM portion

  int nqm;                   // # of QM atoms
  tagint *qmIDs;             // IDs of QM atoms in ascending order
  double **xqm,**fqm;        // QM coords and forces
  double *qqm;               // QM charges
  int *tqm;                  // QM atom types
  double *qpotential;        // Coulomb potential
  double **xqm_mine;         // same values for QM atoms I own
  double *qqm_mine;
  int *tqm_mine;
  double *qpotential_mine;
  int *qm2owned;             // index of local atom for each QM atom
                             // index = -1 if this proc does not own
  
  double *ecoul;             // peratom Coulombic energy from LAMMPS
  int ncoulmax;              // length of ecoul

  // data for MM portion

  int nmm;                  // # of MM atoms
  double **xmm;             // MM coords
  double *qmm;              // MM charges

  // unit conversion factors

  double lmp2mdi_length, mdi2lmp_length;
  double lmp2mdi_energy, mdi2lmp_energy;
  double lmp2mdi_force, mdi2lmp_force;
  double lmp2mdi_pressure, mdi2lmp_pressure;

  // buffers for MDI comm

  int maxbuf;
  int *ibuf1, *ibuf1all;
  double *buf3, *buf3all;

  // methods

  void set_qm2owned();
  void set_qqm();
  void set_tqm();
  void set_xqm();

  void reallocate();
  void send_types();
  void send_elements();
  void send_box();
  void unit_conversions();

  void post_force_direct(int);
  void post_force_potential(int);
  void post_force_aimd(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
