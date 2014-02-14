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

#ifdef COMMAND_CLASS

CommandStyle(prd,PRD)

#else

#ifndef LMP_PRD_H
#define LMP_PRD_H

#include "pointers.h"

namespace LAMMPS_NS {

class PRD : protected Pointers {
 public:
  PRD(class LAMMPS *);
  ~PRD() {}
  void command(int, char **);

 private:
  int me,nprocs;
  int t_event,n_dephase,t_dephase,t_corr;
  double etol,ftol,temp_dephase;
  int maxiter,maxeval,temp_flag,stepmode;
  char *loop_setting,*dist_setting;

  int equal_size_replicas,natoms;
  int neigh_every,neigh_delay,neigh_dist_check;
  int quench_reneighbor;
  bigint nbuild,ndanger;

  double time_dephase,time_dynamics,time_quench,time_comm,time_output;
  double time_start;

  MPI_Comm comm_replica;
  tagint *tagall;
  int *displacements,*imageall;
  double **xall;

  int ncoincident;

  class RanPark *random_select,*random_clock;
  class RanMars *random_dephase;
  class Compute *compute_event;
  class FixEventPRD *fix_event;
  class Velocity *velocity;
  class Compute *temperature;
  class Finish *finish;

  void dephase();
  void dynamics(int, double &);
  void quench();
  int check_event(int replica = -1);
  void share_event(int, int, int);
  void log_event();
  void replicate(int);
  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: PRD command before simulation box is defined

The prd command cannot be used before a read_data,
read_restart, or create_box command.

E: Cannot use PRD with multi-processor replicas unless atom map exists

Use the atom_modify command to create an atom map.

W: Running PRD with only one replica

This is allowed, but you will get no parallel speed-up.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid t_event in prd command

Self-explanatory.

E: PRD nsteps must be multiple of t_event

Self-explanatory.

E: PRD t_corr must be multiple of t_event

Self-explanatory.

E: Could not find compute ID for PRD

Self-explanatory.

W: Resetting reneighboring criteria during PRD

A PRD simulation requires that neigh_modify settings be delay = 0,
every = 1, check = yes.  Since these settings were not in place,
LAMMPS changed them and will restore them to their original values
after the PRD simulation.

E: Too many timesteps

The cummulative timesteps must fit in a 64-bit integer.

E: Cannot use PRD with a changing box

The current box dimensions are not copied between replicas

E: Cannot use PRD with a time-dependent fix defined

PRD alters the timestep in ways that will mess up these fixes.

E: Cannot use PRD with a time-dependent region defined

PRD alters the timestep in ways that will mess up these regions.

E: Cannot use PRD with atom_modify sort enabled

This is a current restriction of PRD.  You must turn off sorting,
which is enabled by default, via the atom_modify command.

E: Too many iterations

You must use a number of iterations that fit in a 32-bit integer
for minimization.

*/
