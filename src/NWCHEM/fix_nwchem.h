/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nwchem,FixNWChem);
// clang-format on
#else

#ifndef LMP_FIX_NWCHEM_H
#define LMP_FIX_NWCHEM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNWChem : public Fix {
 public:
  FixNWChem(class LAMMPS *, int, char **);
  virtual ~FixNWChem();
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
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double compute_scalar() override;
  double memory_usage() override;

 protected:
  char *nw_template;   // template file for NWChem input params
  char *nw_input;      // generated NWChem input file w/ box and coords

  int pbcflag;         // 1 if fully periodic, 0 if fully non-periodic
  int mode;            // AIMD or QMMM
  int qflag;           // 1 if per-atom charge defined, 0 if not
  int qm_init;         // 1 if NWChem and qqm are initialized, 0 if not
  char **elements;     // species name of each LAMMPS atom type

  double qmenergy;     // QM energy
       
  // data for QMMM mode

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

  // conversion factors between LAMMPS and NWChem units

  double lmp2qm_distance,lmp2qm_energy,qm2lmp_force,qm2lmp_energy;

  class Compute *c_pe;      // NOTE: not sure if need this
  class Pair *pair_coul;    // ptr to instance of pair coul/long

  // local methods

  void nwchem_initialize();
  void pre_force_qmmm(int);
  void post_force_qmmm(int);
  void post_force_aimd(int);

  void set_qm2owned();
  void set_qqm();
  void set_tqm();
  void set_xqm();

  void dummy_pspw_input(MPI_Comm, std::string &) {}
  int dummy_pspw_qmmm_minimizer(MPI_Comm, double *, double *, 
                                double *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Must use units metal with fix nwchem command

UNDOCUMENTED

E: Fix nwchem currently runs only in serial

UNDOCUMENTED

E: LAMMPS is linked against incompatible NWChem library

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nwchem does not yet support a LAMMPS calculation of a Coulomb potential

UNDOCUMENTED

E: Could not find fix nwchem compute ID

UNDOCUMENTED

E: Fix nwchem compute ID does not compute pe/atom

UNDOCUMENTED

E: Fix nwchem requires 3d problem

UNDOCUMENTED

E: Fix nwchem cannot compute Coulomb potential

UNDOCUMENTED

E: Fix nwchem requires 3d simulation

UNDOCUMENTED

W: Fix nwchem should come after all other integration fixes

UNDOCUMENTED

E: Internal NWChem problem

UNDOCUMENTED

*/
