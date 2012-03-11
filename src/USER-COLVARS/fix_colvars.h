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
#include <vector>

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
  virtual void min_post_force(int); 
  virtual void post_force(int); 
  virtual void post_force_respa(int, int, int);
  virtual void post_run();
  virtual void end_of_step();
  virtual double compute_scalar();
  virtual double memory_usage();

 protected:
  void   deallocate();        // free internal buffers

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
  int   *taglist;      // list of all atom IDs referenced by colvars.
  std::vector<struct commdata> *coords; // coordinates of colvar atoms
  std::vector<struct commdata> *forces; // received forces of colvar atoms
  std::vector<struct commdata> *oforce; // old total forces of colvar atoms
  
  int   nmax;          // size of atom communication buffer.
  int   size_one;      // bytes per atom in communication buffer.
  struct commdata *comm_buf; // communication buffer
  double *force_buf;   // communication buffer

  void *idmap;         // hash for mapping atom indices to consistent order.
  int  *rev_idmap;     // list of the hash keys for reverse mapping.

  int nlevels_respa;   // flag to determine respa levels.
  int store_forces;    // flag to determine whether to store total forces 
  static  int instances; // count fix instances, since colvars currently
                         // only supports one instance at a time
};

}

#endif
#endif

// Local Variables:
// mode: c++
// compile-command: "make -j4 openmpi"
// c-basic-offset: 2
// fill-column: 76
// indent-tabs-mode: nil
// End:
