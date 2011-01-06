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

#ifdef PAIR_CLASS

PairStyle(reax/c,PairReaxC)

#else

#ifndef LMP_PAIR_REAXC_H
#define LMP_PAIR_REAXC_H

#include "pair.h"
#include "reaxc_types.h"

namespace LAMMPS_NS {

class PairReaxC : public Pair {
 public:
  PairReaxC(class LAMMPS *);
  ~PairReaxC();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void *extract(char *, int &);

 private:
  reax_system *system;
  control_params *control;
  simulation_data *data;
  storage *workspace;
  reax_list *lists;
  output_controls *out_control;
  mpi_datatypes *mpi_data;
  
  double cutmax;
  int *map;
  class FixReaxC *fix_reax;
  
  double *chi,*eta,*gamma;
  int qeqflag;
  int setup_flag;
  int firstwarn;

  void allocate();
  void write_reax_atoms();
  void get_distance(rvec, rvec, double *, rvec *);
  void set_far_nbr(far_neighbor_data *, int, double, rvec);
  int estimate_reax_lists();
  int write_reax_lists();
  void read_reax_forces();
  void setup();
};
  
}

#endif
#endif
