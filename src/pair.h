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
 public:
  double eng_vdwl,eng_coul;      // accumulated energies
  double virial[6];              // accumulated virial

  double cutforce;               // max cutoff for all atom pairs
  double **cutsq;                // max cutoff sq for each atom pair
  int **setflag;                 // 0/1 = whether each i,j has been set

  int comm_forward;              // size of forward communication (0 if none)
  int comm_reverse;              // size of reverse communication (0 if none)

  int single_enable;             // 1 if single() routine exists
  int respa_enable;              // 1 if inner/middle/outer rRESPA routines
  int one_coeff;                 // 1 if allows only one coeff * * call

  int tail_flag;                 // pair_modify flag for LJ tail correction
  double etail,ptail;            // energy/pressure tail corrections
  double etail_ij,ptail_ij;

  class NeighList *list;         // standard neighbor list used by most pairs
  class NeighList *listhalf;     // half list used by some pairs
  class NeighList *listfull;     // full list used by some pairs
  class NeighList *listgranhistory;  // granular history list used by some
  class NeighList *listinner;    // rRESPA lists used by some pairs
  class NeighList *listmiddle;
  class NeighList *listouter;

  struct One {                   // single interaction between 2 atoms
    double fforce;
    double eng_vdwl,eng_coul;
  };

  Pair(class LAMMPS *);
  virtual ~Pair() {}

  // top-level Pair methods

  void init();
  double mix_energy(double, double, double, double);
  double mix_distance(double, double);
  void virial_compute();
  void write_file(int, char **);
  void init_bitmap(double, double, int, int &, int &, int &, int &);
  virtual void modify_params(int, char **);

  // general child-class methods

  virtual void compute(int, int) = 0;
  virtual void compute_inner() {}
  virtual void compute_middle() {}
  virtual void compute_outer(int, int) {}

  virtual void single(int, int, int, int,
		      double, double, double, int, One &) {}

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
  virtual int memory_usage() {return 0;}

  // specific child-class methods for certain Pair styles
  
  virtual void single_embed(int, int, double &) {}

  virtual void *extract_ptr(char *) {return NULL;}

  virtual void extract_charmm(double ***, double ***,
			      double ***, double ***, int *) {}
  virtual void extract_dipole(double ***) {}
  virtual void extract_eam(double *, double **) {}
  virtual void extract_gran(double *, double *, double *, int *) {}
  virtual void extract_long(double *) {}
  virtual void extract_tip4p(double *, int *, int *, int *, int *) {}

 protected:
  int allocated;                       // 0/1 = whether arrays are allocated

                                       // pair_modify settings
  int offset_flag,mix_flag;            // flags for offset and mixing
  int ncoultablebits;                  // size of Coulomb table
  double tabinner;                     // inner cutoff for Coulomb table
};

}

#endif
