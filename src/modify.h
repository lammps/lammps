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

#ifndef MODIFY_H
#define MODIFY_H

#include "stdio.h"
#include "lammps.h"

class Fix;

class Modify : public LAMMPS {
 public:
  int nfix;
  int maxfix;
  int n_initial_integrate,n_pre_decide,n_pre_exchange,n_pre_neighbor;
  int n_post_force,n_final_integrate,n_end_of_step,n_thermo;
  int n_initial_integrate_respa,n_post_force_respa,n_final_integrate_respa;
  int n_min_post_force;
  int nfix_restart_peratom;

  Fix **fix;                 // list of fixes
  int *fmask;                // bit mask for when each fix is applied

  Modify();
  ~Modify();
  void init();
  void setup();
  void initial_integrate();
  void pre_decide();
  void pre_exchange();
  void pre_neighbor();
  void post_force(int);
  void final_integrate();
  void end_of_step();
  int thermo_fields(int, int *, char **);
  void thermo_compute(double *);

  void initial_integrate_respa(int,int);
  void post_force_respa(int,int,int);
  void final_integrate_respa(int);

  void min_post_force(int);

  void add_fix(int, char **);
  void modify_fix(int, char **);
  void delete_fix(char *);

  void write_restart(FILE *);
  int read_restart(FILE *);
  void restart_deallocate();

  int memory_usage();

 private:
                             // lists of fixes to apply at different times
  int *list_initial_integrate,*list_pre_decide;
  int *list_pre_exchange,*list_pre_neighbor;
  int *list_post_force,*list_final_integrate,*list_end_of_step,*list_thermo;
  int *list_initial_integrate_respa,*list_post_force_respa;
  int *list_final_integrate_respa;
  int *list_min_post_force;

  int *end_of_step_every;

  int nfix_restart_global;
  char **id_restart_global;
  char **style_restart_global;
  char **state_restart_global;

  char **id_restart_peratom;
  char **style_restart_peratom;
  int *index_restart_peratom;

  void list_init(int, int &, int *&);
  void list_init_end_of_step(int, int &, int *&);
  void list_init_thermo(int, int &, int *&);
};

#endif
