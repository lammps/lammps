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

#ifndef PRD_H
#define PRD_H

#include "pointers.h"

namespace LAMMPS_NS {

class PRD : protected Pointers {
 public:
  PRD(class LAMMPS *);
  ~PRD();
  void command(int, char **);

 private:
  int me,nprocs;
  int nsteps,t_event,n_dephase,t_dephase,t_corr;
  double etol,ftol,temp_dephase;
  int maxiter,maxeval,temp_flag;
  char *loop_setting,*dist_setting;

  int equal_size_replicas,natoms;
  int neigh_every,neigh_delay,neigh_dist_check;
  int nbuild,ndanger;
  int quench_reneighbor;

  double time_dephase,time_dynamics,time_quench,time_comm,time_output;

  MPI_Comm comm_replica;
  int *tagall,*displacements,*imageall;
  double **xall;
  
  class RanPark *random_select;
  class RanMars *random_dephase;
  class Compute *compute_event;
  class FixEvent *fix_event;
  class Velocity *velocity;
  class Compute *temperature;
  class Finish *finish;

  void dephase();
  void dynamics();
  void quench();
  int check_event(int replica = -1);
  void share_event(int, int);
  void log_event();
  void replicate(int);
  void options(int, char **);
};

}

#endif
