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

#ifndef PAIR_H
#define PAIR_H

#include "pointers.h"

namespace LAMMPS_NS {

class Pair : protected Pointers {
  friend class BondQuartic;
  friend class DihedralCharmm;

 public:
  double eng_vdwl,eng_coul;      // accumulated energies
  double virial[6];              // accumulated virial
  double *eatom,**vatom;         // accumulated per-atom energy/virial

  double cutforce;               // max cutoff for all atom pairs
  double **cutsq;                // max cutoff sq for each atom pair
  int **setflag;                 // 0/1 = whether each i,j has been set

  int comm_forward;              // size of forward communication (0 if none)
  int comm_reverse;              // size of reverse communication (0 if none)

  int single_enable;             // 1 if single() routine exists
  int respa_enable;              // 1 if inner/middle/outer rRESPA routines
  int one_coeff;                 // 1 if allows only one coeff * * call
  int no_virial_compute;         // 1 if does not invoke virial_compute()

  int tail_flag;                 // pair_modify flag for LJ tail correction
  double etail,ptail;            // energy/pressure tail corrections
  double etail_ij,ptail_ij;

  class NeighList *list;         // standard neighbor list used by most pairs
  class NeighList *listhalf;     // half list used by some pairs
  class NeighList *listfull;     // full list used by some pairs
  class NeighList *listgranhistory;  // granular history list used by some pairs
  class NeighList *listinner;    // rRESPA lists used by some pairs
  class NeighList *listmiddle;
  class NeighList *listouter;

  Pair(class LAMMPS *);
  virtual ~Pair();

  // top-level Pair methods

  void init();
  double mix_energy(double, double, double, double);
  double mix_distance(double, double);
  void write_file(int, char **);
  void init_bitmap(double, double, int, int &, int &, int &, int &);
  virtual void modify_params(int, char **);

  // general child-class methods

  virtual void compute(int, int) = 0;
  virtual void compute_inner() {}
  virtual void compute_middle() {}
  virtual void compute_outer(int, int) {}

  virtual double single(int, int, int, int,
			double, double, double, double &) {return 0.0;}

  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;

  virtual void init_style();
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int) {return 0.0;}

  virtual void write_restart(FILE *) {}
  virtual void read_restart(FILE *) {}
  virtual void write_restart_settings(FILE *) {}
  virtual void read_restart_settings(FILE *) {}

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual double memory_usage();

  // specific child-class methods for certain Pair styles
  
  virtual void *extract(char *) {return NULL;}
  virtual void swap_eam(double *, double **) {}
  virtual void reset_dt() {}
  virtual void min_pointers(double **, double **) {}

 protected:
  int allocated;                       // 0/1 = whether arrays are allocated

                                       // pair_modify settings
  int offset_flag,mix_flag;            // flags for offset and mixing
  int ncoultablebits;                  // size of Coulomb table
  double tabinner;                     // inner cutoff for Coulomb table

  // custom data type for accessing Coulomb tables

  typedef union {int i; float f;} union_int_float_t;

  double THIRD;

  int evflag;                          // energy,virial settings
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;
  int vflag_fdotr;
  int maxeatom,maxvatom;

  void ev_setup(int, int);
  void ev_tally(int, int, int, int, double, double, double,
		double, double, double);
  void ev_tally_xyz(int, int, int, int, double, double,
		    double, double, double, double, double, double);
  void ev_tally3(int, int, int, double, double,
		 double *, double *, double *, double *);
  void ev_tally4(int, int, int, int, double,
		 double *, double *, double *, double *, double *, double *);
  void ev_tally_list(int, int *, double, double *);
  void v_tally2(int, int, double, double *);
  void v_tally3(int, int, int, double *, double *, double *, double *);
  void v_tally4(int, int, int, int, double *, double *, double *,
		double *, double *, double *);
  void v_tally_tensor(int, int, int, int,
		      double, double, double, double, double, double);
  void virial_compute();
};

}

#endif
