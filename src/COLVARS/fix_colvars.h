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
      Contributing author:  Axel Kohlmeyer (Temple U)
      Currently maintained by:  Giacomo Fiorin (NIH)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(colvars,FixColvars);
// clang-format on
#else

#ifndef LMP_FIX_COLVARS_H
#define LMP_FIX_COLVARS_H

#include "fix.h"

// Forward declarations
namespace IntHash_NS {
struct inthash_t;
}
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
  char *tfix_name;              // name of thermostat fix.
  int rng_seed;                 // seed to initialize random number generator
  double t_target = 0.0;        // thermostat target temperature
  double energy;                // biasing energy of the fix

  int num_coords;     // total number of atoms controlled by this fix
  tagint *taglist;    // list of all atom IDs referenced by colvars.

  int nmax;                     // size of atom communication buffer.
  int size_one;                 // bytes per atom in communication buffer.
  struct commdata *comm_buf;    // communication buffer
  double *force_buf;            // communication buffer

  /// Arguments passed from fix_modify to the Colvars script interface
  unsigned char *script_args[100];

  IntHash_NS::inthash_t *idmap;    // hash for mapping atom indices to consistent order.

  int nlevels_respa;       // flag to determine respa levels.
  int store_forces;        // flag to determine whether to store total forces
  int unwrap_flag;         // 1 if atom coords are unwrapped, 0 if not
  int init_flag;           // 1 if initialized, 0 if not
  static int instances;    // count fix instances, since colvars currently
                           // only supports one instance at a time
  MPI_Comm root2root;      // inter-root communicator for multi-replica support

  void init_taglist();    // initialize list of atom tags and hash table

  /// Share with Colvars the thermostat fix named by tfix_name
  void set_thermostat_temperature();

  /// Tell Colvars where to get its state from and where to save it
  void setup_io();

  /// Parse LAMMPS-specific arguments to either fix or fix_modify
  /// \param narg Number of arguments
  /// \param arg Array of strings
  /// \param fix_constructor If false, try Colvars commands if LAMMPS ones fail
  int parse_fix_arguments(int narg, char **arg, bool fix_constructor = true);
};

}    // namespace LAMMPS_NS

#endif
#endif
