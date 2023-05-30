// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

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

/* ----------------------------------------------------------------------
      Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(colvars,FixColvars);
// clang-format on
#else

#ifndef LMP_FIX_COLVARS_H
#define LMP_FIX_COLVARS_H

#include "fix.h"

class colvarproxy_lammps;

namespace LAMMPS_NS {

class FixColvars : public Fix {

 public:
  FixColvars(class LAMMPS *, int, char **);
  ~FixColvars() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  int modify_param(int, char **) override;
  void min_setup(int vflag) override { setup(vflag); };
  void min_post_force(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void end_of_step() override;
  void post_run() override;
  double compute_scalar() override;
  double memory_usage() override;

  void write_restart(FILE *) override;
  void restart(char *) override;

 protected:
  colvarproxy_lammps *proxy;    // pointer to the colvars proxy class
  char *conf_file;              // name of colvars config file
  char *inp_name;               // name/prefix of colvars restart file
  char *out_name;               // prefix string for all output files
  char *tmp_name;               // name of thermostat fix.
  int rng_seed;                 // seed to initialize random number generator
  Fix *tstat_fix;               // pointer to thermostat fix
  double energy;                // biasing energy of the fix

  int me;             // my MPI rank in this "world".
  int num_coords;     // total number of atoms controlled by this fix
  tagint *taglist;    // list of all atom IDs referenced by colvars.

  int nmax;                     // size of atom communication buffer.
  int size_one;                 // bytes per atom in communication buffer.
  struct commdata *comm_buf;    // communication buffer
  double *force_buf;            // communication buffer

  void *idmap;       // hash for mapping atom indices to consistent order.
  int *rev_idmap;    // list of the hash keys for reverse mapping.

  int nlevels_respa;       // flag to determine respa levels.
  int store_forces;        // flag to determine whether to store total forces
  int unwrap_flag;         // 1 if atom coords are unwrapped, 0 if not
  int init_flag;           // 1 if initialized, 0 if not
  static int instances;    // count fix instances, since colvars currently
                           // only supports one instance at a time
  MPI_Comm root2root;      // inter-root communicator for multi-replica support
  void one_time_init();    // one time initialization
};

}    // namespace LAMMPS_NS

#endif
#endif
