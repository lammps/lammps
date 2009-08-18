/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef COMPUTE_H
#define COMPUTE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Compute : protected Pointers {
 public:
  char *id,*style;
  int igroup,groupbit;

  double scalar;            // computed global scalar
  double *vector;           // computed global vector
  double *scalar_atom;      // computed per-atom scalar
  double **vector_atom;     // computed per-atom vector

  int scalar_flag;          // 0/1 if compute_scalar() function exists
  int vector_flag;          // 0/1 if compute_vector() function exists
  int size_vector;          // N = size of global vector
  int peratom_flag;         // 0/1 if compute_peratom() function exists
  int size_peratom;         // 0 = scalar_atom, N = size of vector_atom
  int extscalar;            // 0/1 if scalar is intensive/extensive
  int extvector;            // 0/1/-1 if vector is all int/ext/extlist
  int *extlist;             // list of 0/1 int/ext for each vec component

  int tempflag;       // 1 if Compute can be used as temperature
                      // must have both compute_scalar, compute_vector
  int pressflag;      // 1 if Compute can be used as pressure (uses virial)
                      // must have both compute_scalar, compute_vector
  int pressatomflag;  // 1 if Compute calculates per-atom virial
  int peflag;         // 1 if Compute calculates PE (uses Force energies)
  int peatomflag;     // 1 if Compute calculates per-atom PE

  int tempbias;       // 0/1 if Compute temp includes self/extra bias

  int timeflag;       // 1 if Compute stores list of timesteps it's called on
  int ntime;          // # of entries in time list
  int maxtime;        // max # of entries time list can hold
  int *tlist;         // time list of steps the Compute is called on

  int invoked_flag;     // non-zero if invoked or accessed this step, 0 if not
  int invoked_scalar;   // last timestep on which compute_scalar() was invoked
  int invoked_vector;   // ditto for compute_vector()
  int invoked_peratom;  // ditto for compute_peratom()

  double dof;         // degrees-of-freedom for temperature

  int comm_forward;   // size of forward communication (0 if none)
  int comm_reverse;   // size of reverse communication (0 if none)

  Compute(class LAMMPS *, int, char **);
  virtual ~Compute();
  void modify_params(int, char **);
  void reset_extra_dof();

  virtual void init() = 0;
  virtual void init_list(int, class NeighList *) {}
  virtual double compute_scalar() {return 0.0;}
  virtual void compute_vector() {}
  virtual void compute_peratom() {}

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual int dof_remove(int) {return 0;}
  virtual void remove_bias(int, double *) {}
  virtual void remove_bias_all() {}
  virtual void restore_bias(int, double *) {}
  virtual void restore_bias_all() {}

  virtual void reset_extra_compute_fix(char *);

  void addstep(int);
  int matchstep(int);
  void clearstep();

  virtual void *extract(char *) {return NULL;}
  virtual double memory_usage() {return 0.0;}

 protected:
  int extra_dof;               // extra DOF for temperature computes
  int dynamic;                 // recount atoms for temperature computes
  int thermoflag;              // 1 if include fix PE for PE computes

  double vbias[3];             // stored velocity bias for one atom
  double **vbiasall;           // stored velocity bias for all atoms
  int maxbias;                 // size of vbiasall array
};

}

#endif
