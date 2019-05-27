// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

FixStyle(colvars,FixColvars)

#else

#ifndef LMP_FIX_COLVARS_H
#define LMP_FIX_COLVARS_H

#include "fix.h"

// forward declaration
class colvarproxy_lammps;
struct commdata;

namespace LAMMPS_NS {

class FixColvars : public Fix {

 public:
  FixColvars(class LAMMPS *, int, char **);
  virtual ~FixColvars();

  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual int modify_param(int, char **);
  virtual void min_setup(int vflag) {setup(vflag);};
  virtual void min_post_force(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  virtual void end_of_step();
  virtual void post_run();
  virtual double compute_scalar();
  virtual double memory_usage();

  virtual void write_restart(FILE *);
  virtual void restart(char *);

 protected:
  class colvarproxy_lammps *proxy; // pointer to the colvars proxy class
  char *conf_file;     // name of colvars config file
  char *inp_name;      // name/prefix of colvars restart file
  char *out_name;      // prefix string for all output files
  char *tmp_name;      // name of thermostat fix.
  int   rng_seed;      // seed to initialize random number generator
  int   tstat_id;      // id of the thermostat fix
  double energy;       // biasing energy of the fix

  int   me;            // my MPI rank in this "world".
  int   num_coords;    // total number of atoms controlled by this fix
  tagint *taglist;     // list of all atom IDs referenced by colvars.

  int   nmax;          // size of atom communication buffer.
  int   size_one;      // bytes per atom in communication buffer.
  struct commdata *comm_buf; // communication buffer
  double *force_buf;   // communication buffer

  void *idmap;         // hash for mapping atom indices to consistent order.
  int  *rev_idmap;     // list of the hash keys for reverse mapping.

  int nlevels_respa;   // flag to determine respa levels.
  int store_forces;    // flag to determine whether to store total forces
  int unwrap_flag;     // 1 if atom coords are unwrapped, 0 if not
  int init_flag;       // 1 if initialized, 0 if not
  static  int instances; // count fix instances, since colvars currently
                         // only supports one instance at a time
  MPI_Comm root2root;   // inter-root communicator for multi-replica support
  void one_time_init(); // one time initialization
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix colvars for atoms with rmass attribute

The colvars library assigns atom masses per atom type, thus atom styles
which allow setting individual per atom masses are not supported.

E: Missing argument to keyword

Self-explanatory. A keyword was recognized, but no corresponding value
found. Check the input script syntax and compare to the documentation
for the command.

E: Incorrect fix colvars unwrap flag

Self-explanatory. Check the input script syntax.

E: Unknown fix colvars parameter

Self-explanatory. Check your input script syntax.

E: Cannot use fix colvars without atom IDs

Atom IDs are not defined, but fix colvars needs them to identify an atom.

E: Fix colvars requires an atom map, see atom_modify

Use the atom_modify command to create an atom map.

W: Using fix colvars with minimization

Some of the functionality supported with the colvars library is not
meaningful with minimization calculations.

E: Could not find tstat fix ID

Self-explanatory. The thermostat fix ID provided with the tstat keyword
is not defined (yet or anymore). Check your input file.

E: Run aborted on request from colvars module

Some error condition happened inside the colvars library that prohibits
it from continuing. Please examine the output for additional information.

*/

