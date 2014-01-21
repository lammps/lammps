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

/* ----------------------------------------------------------------------
   Contributing authors: Mark Stevens (SNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "respa.h"
#include "neighbor.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix_respa.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Respa::Respa(LAMMPS *lmp, int narg, char **arg) : Integrate(lmp, narg, arg)
{
  if (narg < 1) error->all(FLERR,"Illegal run_style respa command");

  nlevels = force->inumeric(FLERR,arg[0]);
  if (nlevels < 1) error->all(FLERR,"Respa levels must be >= 1");

  if (narg < nlevels) error->all(FLERR,"Illegal run_style respa command");
  loop = new int[nlevels];
  for (int iarg = 1; iarg < nlevels; iarg++) {
    loop[iarg-1] = force->inumeric(FLERR,arg[iarg]);
    if (loop[iarg-1] <= 0) error->all(FLERR,"Illegal run_style respa command");
  }
  loop[nlevels-1] = 1;

  // set level at which each force is computed
  // argument settings override defaults

  level_bond = level_angle = level_dihedral = level_improper = -1;
  level_pair = level_kspace = -1;
  level_inner = level_middle = level_outer = -1;

  int iarg = nlevels;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_bond = force->inumeric(FLERR,arg[iarg+1]) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_angle = force->inumeric(FLERR,arg[iarg+1]) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_dihedral = force->inumeric(FLERR,arg[iarg+1]) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"improper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_improper = force->inumeric(FLERR,arg[iarg+1]) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_pair = force->inumeric(FLERR,arg[iarg+1]) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"inner") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_inner = force->inumeric(FLERR,arg[iarg+1]) - 1;
      cutoff[0] = force->numeric(FLERR,arg[iarg+2]);
      cutoff[1] = force->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"middle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_middle = force->inumeric(FLERR,arg[iarg+1]) - 1;
      cutoff[2] = force->numeric(FLERR,arg[iarg+2]);
      cutoff[3] = force->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"outer") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_outer = force->inumeric(FLERR,arg[iarg+1]) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_kspace = force->inumeric(FLERR,arg[iarg+1]) - 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal run_style respa command");
  }

  // cannot specify both pair and inner/middle/outer

  if (level_pair >= 0 &&
      (level_inner >= 0 || level_middle >= 0 || level_outer >= 0))
    error->all(FLERR,"Cannot set both respa pair and inner/middle/outer");

  // if either inner and outer is specified, then both must be

  if ((level_inner >= 0 && level_outer == -1) ||
      (level_outer >= 0 && level_inner == -1))
    error->all(FLERR,"Must set both respa inner and outer");

  // middle cannot be set without inner/outer

  if (level_middle >= 0 && level_inner == -1)
    error->all(FLERR,"Cannot set respa middle without inner/outer");

  // set defaults if user did not specify level
  // bond to innermost level
  // angle same as bond, dihedral same as angle, improper same as dihedral
  // pair to outermost level if no inner/middle/outer
  // inner/middle/outer have no defaults
  // kspace same as pair or outer

  if (level_bond == -1) level_bond = 0;
  if (level_angle == -1) level_angle = level_bond;
  if (level_dihedral == -1) level_dihedral = level_angle;
  if (level_improper == -1) level_improper = level_dihedral;
  if (level_pair == -1 && level_inner == -1) level_pair = nlevels-1;
  if (level_kspace == -1 && level_pair >= 0) level_kspace = level_pair;
  if (level_kspace == -1 && level_pair == -1) level_kspace = level_outer;

  // print respa levels

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Respa levels:\n");
      for (int i = 0; i < nlevels; i++) {
        fprintf(screen,"  %d =",i+1);
        if (level_bond == i) fprintf(screen," bond");
        if (level_angle == i) fprintf(screen," angle");
        if (level_dihedral == i) fprintf(screen," dihedral");
        if (level_improper == i) fprintf(screen," improper");
        if (level_pair == i) fprintf(screen," pair");
        if (level_inner == i) fprintf(screen," pair-inner");
        if (level_middle == i) fprintf(screen," pair-middle");
        if (level_outer == i) fprintf(screen," pair-outer");
        if (level_kspace == i) fprintf(screen," kspace");
        fprintf(screen,"\n");
      }
    }
    if (logfile) {
      fprintf(logfile,"Respa levels:\n");
      for (int i = 0; i < nlevels; i++) {
        fprintf(logfile,"  %d =",i+1);
        if (level_bond == i) fprintf(logfile," bond");
        if (level_angle == i) fprintf(logfile," angle");
        if (level_dihedral == i) fprintf(logfile," dihedral");
        if (level_improper == i) fprintf(logfile," improper");
        if (level_pair == i) fprintf(logfile," pair");
        if (level_inner == i) fprintf(logfile," pair-inner");
        if (level_middle == i) fprintf(logfile," pair-middle");
        if (level_outer == i) fprintf(logfile," pair-outer");
        if (level_kspace == i) fprintf(logfile," kspace");
        fprintf(logfile,"\n");
      }
    }
  }

  // check that levels are in correct order

  if (level_angle < level_bond || level_dihedral < level_angle ||
      level_improper < level_dihedral)
    error->all(FLERR,"Invalid order of forces within respa levels");
  if (level_pair >= 0) {
    if (level_pair < level_improper || level_kspace < level_pair)
      error->all(FLERR,"Invalid order of forces within respa levels");
  }
  if (level_pair == -1 && level_middle == -1) {
    if (level_inner < level_improper || level_outer < level_inner ||
        level_kspace < level_outer)
      error->all(FLERR,"Invalid order of forces within respa levels");
  }
  if (level_pair == -1 && level_middle >= 0) {
    if (level_inner < level_improper || level_middle < level_inner ||
        level_outer < level_inner || level_kspace < level_outer)
      error->all(FLERR,"Invalid order of forces within respa levels");
  }

  // warn if any levels are devoid of forces

  int flag = 0;
  for (int i = 0; i < nlevels; i++)
    if (level_bond != i && level_angle != i && level_dihedral != i &&
        level_improper != i && level_pair != i && level_inner != i &&
        level_middle != i && level_outer != i && level_kspace != i) flag = 1;
  if (flag && comm->me == 0)
    error->warning(FLERR,"One or more respa levels compute no forces");

  // check cutoff consistency if inner/middle/outer are enabled

  if (level_inner >= 0 && cutoff[1] < cutoff[0])
    error->all(FLERR,"Respa inner cutoffs are invalid");
  if (level_middle >= 0 && (cutoff[3] < cutoff[2] || cutoff[2] < cutoff[1]))
    error->all(FLERR,"Respa middle cutoffs are invalid");

  // set outer pair of cutoffs to inner pair if middle is not enabled

  if (level_inner >= 0 && level_middle < 0) {
    cutoff[2] = cutoff[0];
    cutoff[3] = cutoff[1];
  }

  // allocate other needed arrays

  newton = new int[nlevels];
  step = new double[nlevels];
}

/* ---------------------------------------------------------------------- */

Respa::~Respa()
{
  delete [] loop;
  delete [] newton;
  delete [] step;
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void Respa::init()
{
  Integrate::init();

  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning(FLERR,"No fixes defined, atoms won't move");

  // create fix needed for storing atom-based respa level forces
  // will delete it at end of run

  char **fixarg = new char*[4];
  fixarg[0] = (char *) "RESPA";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "RESPA";
  fixarg[3] = new char[8];
  sprintf(fixarg[3],"%d",nlevels);
  modify->add_fix(4,fixarg);
  delete [] fixarg[3];
  delete [] fixarg;
  fix_respa = (FixRespa *) modify->fix[modify->nfix-1];

  // insure respa inner/middle/outer is using Pair class that supports it

  if (level_inner >= 0)
    if (force->pair && force->pair->respa_enable == 0)
      error->all(FLERR,"Pair style does not support rRESPA inner/middle/outer");

  // virial_style = 1 (explicit) since never computed implicitly like Verlet

  virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // detect if fix omp is present and will clear force arrays

  int ifix = modify->find_fix("package_omp");
  if (ifix >= 0) external_force_clear = 1;

  // set flags for what arrays to clear in force_clear()
  // need to clear additionals arrays if they exist

  torqueflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  erforceflag = 0;
  if (atom->erforce_flag) erforceflag = 1;
  e_flag = 0;
  if (atom->e_flag) e_flag = 1;
  rho_flag = 0;
  if (atom->rho_flag) rho_flag = 1;

  // step[] = timestep for each level

  step[nlevels-1] = update->dt;
  for (int ilevel = nlevels-2; ilevel >= 0; ilevel--)
    step[ilevel] = step[ilevel+1]/loop[ilevel];

  // set newton flag for each level

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    newton[ilevel] = 0;
    if (force->newton_bond) {
      if (level_bond == ilevel || level_angle == ilevel ||
          level_dihedral == ilevel || level_improper == ilevel)
        newton[ilevel] = 1;
    }
    if (force->newton_pair) {
      if (level_pair == ilevel || level_inner == ilevel ||
          level_middle == ilevel || level_outer == ilevel)
        newton[ilevel] = 1;
    }
  }

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void Respa::setup()
{
  if (comm->me == 0 && screen) fprintf(screen,"Setting up run ...\n");

  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atom->setup();
  modify->setup_pre_exchange();
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();
  neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear(newton[ilevel]);
    modify->setup_pre_force_respa(vflag,ilevel);
    if (level_pair == ilevel && pair_compute_flag)
      force->pair->compute(eflag,vflag);
    if (level_inner == ilevel && pair_compute_flag)
      force->pair->compute_inner();
    if (level_middle == ilevel && pair_compute_flag)
      force->pair->compute_middle();
    if (level_outer == ilevel && pair_compute_flag)
      force->pair->compute_outer(eflag,vflag);
    if (level_bond == ilevel && force->bond)
      force->bond->compute(eflag,vflag);
    if (level_angle == ilevel && force->angle)
      force->angle->compute(eflag,vflag);
    if (level_dihedral == ilevel && force->dihedral)
      force->dihedral->compute(eflag,vflag);
    if (level_improper == ilevel && force->improper)
      force->improper->compute(eflag,vflag);
    if (level_kspace == ilevel && force->kspace) {
      force->kspace->setup();
      if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    }
    if (newton[ilevel]) comm->reverse_comm();
    copy_f_flevel(ilevel);
  }

  sum_flevel_f();
  modify->setup(vflag);
  output->setup();
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void Respa::setup_minimal(int flag)
{
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    modify->setup_pre_exchange();
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();
    neighbor->build();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear(newton[ilevel]);
    modify->setup_pre_force_respa(vflag,ilevel);
    if (level_pair == ilevel && pair_compute_flag)
      force->pair->compute(eflag,vflag);
    if (level_inner == ilevel && pair_compute_flag)
      force->pair->compute_inner();
    if (level_middle == ilevel && pair_compute_flag)
      force->pair->compute_middle();
    if (level_outer == ilevel && pair_compute_flag)
      force->pair->compute_outer(eflag,vflag);
    if (level_bond == ilevel && force->bond)
      force->bond->compute(eflag,vflag);
    if (level_angle == ilevel && force->angle)
      force->angle->compute(eflag,vflag);
    if (level_dihedral == ilevel && force->dihedral)
      force->dihedral->compute(eflag,vflag);
    if (level_improper == ilevel && force->improper)
      force->improper->compute(eflag,vflag);
    if (level_kspace == ilevel && force->kspace) {
      force->kspace->setup();
      if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    }
    if (newton[ilevel]) comm->reverse_comm();
    copy_f_flevel(ilevel);
  }

  sum_flevel_f();
  modify->setup(vflag);
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void Respa::run(int n)
{
  bigint ntimestep;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    recurse(nlevels-1);
    sum_flevel_f();           // needed in case end_of_step or output use f

    if (modify->n_end_of_step) modify->end_of_step();

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(update->ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}

/* ----------------------------------------------------------------------
   delete rRESPA fix at end of run, so its atom arrays won't persist
------------------------------------------------------------------------- */

void Respa::cleanup()
{
  modify->post_run();
  modify->delete_fix("RESPA");
  domain->box_too_small_check();
  update->update_time();
}

/* ---------------------------------------------------------------------- */

void Respa::reset_dt()
{
  step[nlevels-1] = update->dt;
  for (int ilevel = nlevels-2; ilevel >= 0; ilevel--)
    step[ilevel] = step[ilevel+1]/loop[ilevel];
}

/* ---------------------------------------------------------------------- */

void Respa::recurse(int ilevel)
{
  copy_flevel_f(ilevel);

  for (int iloop = 0; iloop < loop[ilevel]; iloop++) {

    modify->initial_integrate_respa(vflag,ilevel,iloop);
    if (modify->n_post_integrate_respa)
      modify->post_integrate_respa(ilevel,iloop);

    if (ilevel) recurse(ilevel-1);

    // at outermost level, check on rebuilding neighbor list
    // at innermost level, communicate
    // at middle levels, do nothing

    if (ilevel == nlevels-1) {
      int nflag = neighbor->decide();
      if (nflag) {
        if (modify->n_pre_exchange) modify->pre_exchange();
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
          domain->reset_box();
          comm->setup();
          if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (atom->sortfreq > 0 &&
            update->ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(TIME_COMM);
        if (modify->n_pre_neighbor) modify->pre_neighbor();
        neighbor->build();
        timer->stamp(TIME_NEIGHBOR);
      } else if (ilevel == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(TIME_COMM);
      }

    } else if (ilevel == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(TIME_COMM);
    }

    // force computations
    // important that ordering is same as Verlet
    // so that any order dependencies are the same
    // when potentials are invoked at same level

    force_clear(newton[ilevel]);
    if (modify->n_pre_force_respa)
      modify->pre_force_respa(vflag,ilevel,iloop);

    timer->stamp();
    if (level_pair == ilevel && pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }
    if (level_inner == ilevel && pair_compute_flag) {
      force->pair->compute_inner();
      timer->stamp(TIME_PAIR);
    }
    if (level_middle == ilevel && pair_compute_flag) {
      force->pair->compute_middle();
      timer->stamp(TIME_PAIR);
    }
    if (level_outer == ilevel && pair_compute_flag) {
      force->pair->compute_outer(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }
    if (level_bond == ilevel && force->bond) {
      force->bond->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_angle == ilevel && force->angle) {
      force->angle->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_dihedral == ilevel && force->dihedral) {
      force->dihedral->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_improper == ilevel && force->improper) {
      force->improper->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_kspace == ilevel && kspace_compute_flag) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    }

    if (newton[ilevel]) {
      comm->reverse_comm();
      timer->stamp(TIME_COMM);
    }

    if (modify->n_post_force_respa)
      modify->post_force_respa(vflag,ilevel,iloop);
    modify->final_integrate_respa(ilevel,iloop);
  }

  copy_f_flevel(ilevel);
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
------------------------------------------------------------------------- */

void Respa::force_clear(int newtonflag)
{
  if (external_force_clear) return;

  // clear global force array
  // nall includes ghosts only if newton flag is set

  int nall;
  if (newtonflag) nall = atom->nlocal + atom->nghost;
  else nall = atom->nlocal;

  size_t nbytes = sizeof(double) * nall;

  if (nbytes > 0 ) {
    memset(&(atom->f[0][0]),0,3*nbytes);
    if (torqueflag)  memset(&(atom->torque[0][0]),0,3*nbytes);
    if (erforceflag) memset(&(atom->erforce[0]),  0,  nbytes);
    if (e_flag)      memset(&(atom->de[0]),       0,  nbytes);
    if (rho_flag)    memset(&(atom->drho[0]),     0,  nbytes);
  }
}

/* ----------------------------------------------------------------------
   copy force components from atom->f to FixRespa->f_level
------------------------------------------------------------------------- */

void Respa::copy_f_flevel(int ilevel)
{
  double ***f_level = fix_respa->f_level;
  double **f = atom->f;
  int n = atom->nlocal;

  for (int i = 0; i < n; i++) {
    f_level[i][ilevel][0] = f[i][0];
    f_level[i][ilevel][1] = f[i][1];
    f_level[i][ilevel][2] = f[i][2];
  }
}

/* ----------------------------------------------------------------------
   copy force components from FixRespa->f_level to atom->f
------------------------------------------------------------------------- */

void Respa::copy_flevel_f(int ilevel)
{
  double ***f_level = fix_respa->f_level;
  double **f = atom->f;
  int n = atom->nlocal;

  for (int i = 0; i < n; i++) {
    f[i][0] = f_level[i][ilevel][0];
    f[i][1] = f_level[i][ilevel][1];
    f[i][2] = f_level[i][ilevel][2];
  }
}

/* ----------------------------------------------------------------------
   sum all force components from FixRespa->f_level to create full atom->f
------------------------------------------------------------------------- */

void Respa::sum_flevel_f()
{
  copy_flevel_f(0);

  double ***f_level = fix_respa->f_level;
  double **f = atom->f;
  int n = atom->nlocal;

  for (int ilevel = 1; ilevel < nlevels; ilevel++) {
    for (int i = 0; i < n; i++) {
      f[i][0] += f_level[i][ilevel][0];
      f[i][1] += f_level[i][ilevel][1];
      f[i][2] += f_level[i][ilevel][2];
    }
  }
}
