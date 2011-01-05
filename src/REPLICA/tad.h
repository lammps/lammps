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

CommandStyle(tad,TAD)

#else

#ifndef LMP_TAD_H
#define LMP_TAD_H

#include "pointers.h"

namespace LAMMPS_NS {

class TAD : protected Pointers {
 public:
  TAD(class LAMMPS *);
  ~TAD();
  void command(int, char **);

 private:
  int me,nprocs;
  int nsteps,t_event;
  double templo,temphi,delta_conf,tmax;
  double etol,ftol,etol_neb,ftol_neb;
  int maxiter,maxeval,n1steps_neb,n2steps_neb,nevery_neb;
  char *min_style, *min_style_neb;
  double delta_beta,ratio_beta;
  double deltconf,deltstop,deltfirst; // Times since last event
  int event_first;

  int neigh_every,neigh_delay,neigh_dist_check;
  int nbuild,ndanger;
  int quench_reneighbor;

  double time_dynamics,time_quench,time_neb,time_comm,time_output;
  double time_start;

  class NEB *neb;                    // NEB object
  class Fix *fix_neb;                 // FixNEB object
  class Compute *compute_event;      // compute to detect event
  class FixEventTAD *fix_event;      // current event/state
  class FixStoreState *fix_revert;   // revert state
  FixEventTAD **fix_event_list;      // list of possible events
  int n_event_list;                  // number of events
  int nmax_event_list;               // allocated events
  int nmin_event_list;               // minimum allocation

  char *neb_logfilename;             // filename for ulogfile_neb
  FILE *uscreen_neb;                 // neb universe screen output
  FILE *ulogfile_neb;                // neb universe logfile
  FILE *uscreen_lammps;              // lammps universe screen output
  FILE *ulogfile_lammps;             // lammps universe logfile

  class Finish *finish;

  void dynamics();
  void quench();
  int check_event();
  int check_confidence();
  void perform_neb(int);
  void log_event();
  void options(int, char **);
  void revert();
  void add_event();
  void perform_event(int);
  void compute_tlo(int);
  void grow_event_list(int);
  void initialize_event_list();
  void delete_event_list();
};

}

#endif
#endif
