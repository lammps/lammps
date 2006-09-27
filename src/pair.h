/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_H
#define PAIR_H

#include "lammps.h"

class Pair : public LAMMPS {
 public:
                                       // used by other classes
  double cutforce;                     // max cutoff for all atom pairs
  double eng_vdwl,eng_coul;            // accumulated energies
  double virial[6];                    // accumulated virial

                                       // used by neighbor class
  int neigh_half_every;                // 0/1 = if needs half neighbor list
  int neigh_full_every;                // 0/1 = if needs full neighbor list

                                       // common to all pair classes
  int allocated;                       // 0/1 = whether arrays are allocated
  int **setflag;                       // 0/1 = whether each i,j has been set
  double **cutsq;                      // max cutoff sq for each atom pair

                                       // pair modify settings
  int offset_flag,mix_flag;            // flags for offset and mixing
  int tail_flag;                       // flag for LJ tail correction
  double etail,ptail;                  // energy/pressure tail corrections
  double etail_ij,ptail_ij;
  int ncoultablebits;                  // size of Coulomb table
  double tabinner;                     // inner cutoff for Coulomb table

  int single_enable;                   // 0/1 if single() routine exists
  int respa_enable;                    // 0/1 if inner/middle/outer routines
  int one_coeff;                       // 0/1 if allows only one coeff * * call

  struct One {                         // single interaction between 2 atoms
    double fforce;
    double eng_vdwl,eng_coul;
    double fx,fy,fz;
    double tix,tiy,tiz,tjx,tjy,tjz;
  };

  Pair();
  virtual ~Pair() {}

  void init();
  double mix_energy(double, double, double, double);
  double mix_distance(double, double);
  void virial_compute();
  void write_file(int, char **);
  void init_bitmap(double, double, int, int &, int &, int &, int &);

  virtual void modify_params(int, char **);

  virtual void compute(int, int) = 0;
  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;

  virtual double init_one(int, int) {return 0.0;}
  virtual void init_style() {}
  virtual void compute_inner() {}
  virtual void compute_middle() {}
  virtual void compute_outer(int, int) {}
  virtual void single(int, int, int, int,
		      double, double, double, int, One &) {}
  virtual void single_embed(int, int, double &, int, double &) {}

  virtual void write_restart(FILE *) {}
  virtual void read_restart(FILE *) {}
  virtual void write_restart_settings(FILE *) {}
  virtual void read_restart_settings(FILE *) {}

  virtual int pack_comm(int, int *, double *, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual int memory_usage() {return 0;}
};

#endif
