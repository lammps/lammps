// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "respa.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_respa.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "timer.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Respa::Respa(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg),
  step(nullptr), loop(nullptr), hybrid_level(nullptr), hybrid_compute(nullptr),
  newton(nullptr), fix_respa(nullptr)
{
  nhybrid_styles = 0;
  if (narg < 1) error->all(FLERR,"Illegal run_style respa command");

  nlevels = utils::inumeric(FLERR,arg[0],false,lmp);
  if (nlevels < 1) error->all(FLERR,"Respa levels must be >= 1");

  if (narg < nlevels) error->all(FLERR,"Illegal run_style respa command");
  loop = new int[nlevels];
  for (int iarg = 1; iarg < nlevels; iarg++) {
    loop[iarg-1] = utils::inumeric(FLERR,arg[iarg],false,lmp);
    if (loop[iarg-1] <= 0) error->all(FLERR,"Illegal run_style respa command");
  }
  loop[nlevels-1] = 1;

  // set level at which each force is computed
  // argument settings override defaults

  level_bond = level_angle = level_dihedral = level_improper = -1;
  level_pair = level_kspace = -1;
  level_inner = level_middle = level_outer = -1;

  // defaults for hybrid pair styles

  nhybrid_styles = 0;
  tally_global = 1;
  pair_compute = 1;

  int iarg = nlevels;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_bond = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_angle = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_dihedral = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"improper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_improper = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_pair = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"inner") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_inner = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      cutoff[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      cutoff[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"middle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_middle = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      cutoff[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      cutoff[3] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"outer") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_outer = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run_style respa command");
      level_kspace = utils::inumeric(FLERR,arg[iarg+1],false,lmp) - 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"hybrid") == 0) {
      // the hybrid keyword requires a hybrid pair style
      if (!utils::strmatch(force->pair_style,"^hybrid"))
        error->all(FLERR,"Illegal run_style respa command");
      PairHybrid *hybrid = (PairHybrid *) force->pair;
      nhybrid_styles = hybrid->nstyles;
      // each hybrid sub-style needs to be assigned to a respa level
      if (iarg+nhybrid_styles > narg)
        error->all(FLERR,"Illegal run_style respa command");
      hybrid_level = new int[nhybrid_styles];
      hybrid_compute = new int[nhybrid_styles];
      for (int i=0; i < nhybrid_styles; ++i) {
        ++iarg;
        hybrid_level[i] = utils::inumeric(FLERR,arg[iarg],false,lmp)-1;
      }
      ++iarg;
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

  // cannot combine hybrid with any of pair/inner/middle/outer

  if ((nhybrid_styles > 0) && (level_pair >= 0 || level_inner >= 0
                               || level_middle >= 0 || level_outer >= 0))
    error->all(FLERR,"Cannot set respa hybrid and "
               "any of pair/inner/middle/outer");

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

  if (level_pair == -1 && level_inner == -1 && nhybrid_styles < 1)
    level_pair = nlevels-1;

  if (level_kspace == -1 && level_pair >= 0) level_kspace = level_pair;
  if (level_kspace == -1 && level_pair == -1) {
    if (nhybrid_styles < 1) {
      level_kspace = level_outer;
    } else {
      int max_hybrid_level = -1;
      for (int i=0; i < nhybrid_styles; ++i) {
        if (max_hybrid_level < hybrid_level[i])
          max_hybrid_level = hybrid_level[i];
      }
      level_kspace = max_hybrid_level;
    }
  }

  // print respa levels

  if (comm->me == 0) {
    std::string mesg = "Respa levels:\n";
    for (int i = 0; i < nlevels; i++) {
      mesg += fmt::format("  {} =",i+1);
      if (level_bond == i)      mesg += " bond";
      if (level_angle == i)     mesg += " angle";
      if (level_dihedral == i)  mesg += " dihedral";
      if (level_improper == i)  mesg += " improper";
      if (level_pair == i)      mesg += " pair";
      if (level_inner == i)     mesg += " pair-inner";
      if (level_middle == i)    mesg += " pair-middle";
      if (level_outer == i)     mesg += " pair-outer";
      for (int j=0; j < nhybrid_styles; j++)
        if (hybrid_level[j] == i) mesg += fmt::format(" hybrid-{}",j+1);
      if (level_kspace == i)    mesg += " kspace";
      mesg += "\n";
    }
    utils::logmesg(lmp,mesg);
  }

  // check that levels are in correct order

  if (level_angle < level_bond || level_dihedral < level_angle ||
      level_improper < level_dihedral)
    error->all(FLERR,"Invalid order of forces within respa levels");
  if (level_pair >= 0) {
    if (level_pair < level_improper || level_kspace < level_pair)
      error->all(FLERR,"Invalid order of forces within respa levels");
  }
  if (level_pair == -1 && level_middle == -1 && nhybrid_styles < 1) {
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

  // ensure that pair->compute() is run properly
  // when the hybrid keyword is not used

  if (nhybrid_styles < 1) {
    pair_compute = 1;
    tally_global = 1;
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
  if (nhybrid_styles > 0) {
    delete [] hybrid_level;
    delete [] hybrid_compute;
  }
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
  // if supported, we also store torques on a per-level basis

  std::string cmd = fmt::format("RESPA all RESPA {}",nlevels);
  if (atom->torque_flag) cmd += " torque";
  fix_respa = (FixRespa *) modify->add_fix(cmd);

  // insure respa inner/middle/outer is using Pair class that supports it

  if (level_inner >= 0)
    if (force->pair && force->pair->respa_enable == 0)
      error->all(FLERR,"Pair style does not support rRESPA inner/middle/outer");

  // virial_style = VIRIAL_PAIR (explicit)
  //   since never computed implicitly with virial_fdotr_compute() like Verlet

  virial_style = VIRIAL_PAIR;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // detect if fix omp is present and will clear force arrays

  int ifix = modify->find_fix("package_omp");
  if (ifix >= 0) external_force_clear = 1;

  // set flags for arrays to clear in force_clear()

  torqueflag = extraflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  if (atom->avec->forceclearflag) extraflag = 1;

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

      if (nhybrid_styles > 0) {
        set_compute_flags(ilevel);
        if (pair_compute) newton[ilevel] = 1;
      }
    }
  }

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void Respa::setup(int flag)
{
  if (comm->me == 0 && screen) {
    std::string mesg = "Setting up r-RESPA run ...\n";
    if (flag) {
      mesg += fmt::format("  Unit style    : {}\n",update->unit_style);
      mesg += fmt::format("  Current step  : {}\n", update->ntimestep);

      mesg += "  Time steps    :";
      for (int ilevel=0; ilevel < nlevels; ++ilevel)
        mesg += fmt::format(" {}:{}",ilevel+1, step[ilevel]);

      mesg += "\n  r-RESPA fixes :";
      for (int l=0; l < modify->n_post_force_respa; ++l) {
        Fix *f = modify->fix[modify->list_post_force_respa[l]];
        if (f->respa_level >= 0)
          mesg += fmt::format(" {}:{}[{}]",
                              MIN(f->respa_level+1,nlevels),
                              f->style,f->id);
      }
      mesg += "\n";
      fputs(mesg.c_str(),screen);
      timer->print_timeout(screen);
    }
  }

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
  neighbor->build(1);
  modify->setup_post_neighbor();
  neighbor->ncalls = 0;

  // compute all forces

  force->setup();
  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear(newton[ilevel]);
    modify->setup_pre_force_respa(vflag,ilevel);

    if (nhybrid_styles > 0) {
      set_compute_flags(ilevel);
      force->pair->compute(eflag,vflag);
    }
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

    modify->setup_pre_reverse(eflag,vflag);
    if (newton[ilevel]) comm->reverse_comm();
    copy_f_flevel(ilevel);
  }

  sum_flevel_f();
  modify->setup(vflag);
  output->setup(flag);
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
    neighbor->build(1);
    modify->setup_post_neighbor();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear(newton[ilevel]);
    modify->setup_pre_force_respa(vflag,ilevel);

    if (nhybrid_styles > 0) {
      set_compute_flags(ilevel);
      force->pair->compute(eflag,vflag);
    }

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

    modify->setup_pre_reverse(eflag,vflag);
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
    if (timer->check_timeout(i)) {
      update->nsteps = i;
      break;
    }

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    recurse(nlevels-1);

    // needed in case end_of_step() or output() use total force

    sum_flevel_f();

    if (modify->n_end_of_step) {
      timer->stamp();
      modify->end_of_step();
      timer->stamp(Timer::MODIFY);
    }

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(update->ntimestep);
      timer->stamp(Timer::OUTPUT);
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

    timer->stamp();
    modify->initial_integrate_respa(vflag,ilevel,iloop);
    if (modify->n_post_integrate_respa)
      modify->post_integrate_respa(ilevel,iloop);
    timer->stamp(Timer::MODIFY);

    // at outermost level, check on rebuilding neighbor list
    // at innermost level, communicate
    // at middle levels, do nothing

    if (ilevel == nlevels-1) {
      int nflag = neighbor->decide();
      if (nflag) {
        if (modify->n_pre_exchange) {
          timer->stamp();
          modify->pre_exchange();
          timer->stamp(Timer::MODIFY);
        }
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
        timer->stamp(Timer::COMM);
        if (modify->n_pre_neighbor) {
          modify->pre_neighbor();
          timer->stamp(Timer::MODIFY);
        }
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
        if (modify->n_post_neighbor) {
          modify->post_neighbor();
          timer->stamp(Timer::MODIFY);
        }

      } else if (ilevel == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
      }

    } else if (ilevel == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(Timer::COMM);
    }

    // rRESPA recursion thru all levels
    // this used to be before neigh list build,
    // which prevented per-atom energy/stress being tallied correctly
    // b/c atoms migrated to new procs between short/long force calls
    // now they migrate at very start of rRESPA timestep, before all forces

    if (ilevel) recurse(ilevel-1);

    // force computations
    // important that ordering is same as Verlet
    // so that any order dependencies are the same
    // when potentials are invoked at same level

    force_clear(newton[ilevel]);
    if (modify->n_pre_force_respa) {
      timer->stamp();
      modify->pre_force_respa(vflag,ilevel,iloop);
      timer->stamp(Timer::MODIFY);
    }

    timer->stamp();
    if (nhybrid_styles > 0) {
      set_compute_flags(ilevel);
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_pair == ilevel && pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_inner == ilevel && pair_compute_flag) {
      force->pair->compute_inner();
      timer->stamp(Timer::PAIR);
    }
    if (level_middle == ilevel && pair_compute_flag) {
      force->pair->compute_middle();
      timer->stamp(Timer::PAIR);
    }
    if (level_outer == ilevel && pair_compute_flag) {
      force->pair->compute_outer(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_bond == ilevel && force->bond) {
      force->bond->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_angle == ilevel && force->angle) {
      force->angle->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_dihedral == ilevel && force->dihedral) {
      force->dihedral->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_improper == ilevel && force->improper) {
      force->improper->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_kspace == ilevel && kspace_compute_flag) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(Timer::KSPACE);
    }

    if (modify->n_pre_reverse) {
      modify->pre_reverse(eflag,vflag);
      timer->stamp(Timer::MODIFY);
    }

    if (newton[ilevel]) {
      comm->reverse_comm();
      timer->stamp(Timer::COMM);
    }
    timer->stamp();
    if (modify->n_post_force_respa)
      modify->post_force_respa(vflag,ilevel,iloop);
    modify->final_integrate_respa(ilevel,iloop);
    timer->stamp(Timer::MODIFY);
  }

  copy_f_flevel(ilevel);
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void Respa::force_clear(int /*newtonflag*/)
{
  if (external_force_clear) return;

  // clear global force array
  // if either newton flag is set, also include ghosts

  size_t nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  if (nbytes) {
    memset(&atom->f[0][0],0,3*nbytes);
    if (torqueflag) memset(&atom->torque[0][0],0,3*nbytes);
    if (extraflag) atom->avec->force_clear(0,nbytes);
  }
}

/* ----------------------------------------------------------------------
   copy force components from atom->f to FixRespa->f_level
------------------------------------------------------------------------- */

void Respa::copy_f_flevel(int ilevel)
{
  double ***f_level = fix_respa->f_level;
  double **f = atom->f;
  double ***t_level = fix_respa->t_level;
  double **t = atom->torque;
  int n = atom->nlocal;

  for (int i = 0; i < n; i++) {
    f_level[i][ilevel][0] = f[i][0];
    f_level[i][ilevel][1] = f[i][1];
    f_level[i][ilevel][2] = f[i][2];
    if (fix_respa->store_torque) {
      t_level[i][ilevel][0] = t[i][0];
      t_level[i][ilevel][1] = t[i][1];
      t_level[i][ilevel][2] = t[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   copy force components from FixRespa->f_level to atom->f
------------------------------------------------------------------------- */

void Respa::copy_flevel_f(int ilevel)
{
  double ***f_level = fix_respa->f_level;
  double **f = atom->f;
  double ***t_level = fix_respa->t_level;
  double **t = atom->torque;
  int n = atom->nlocal;

  for (int i = 0; i < n; i++) {
    f[i][0] = f_level[i][ilevel][0];
    f[i][1] = f_level[i][ilevel][1];
    f[i][2] = f_level[i][ilevel][2];
    if (fix_respa->store_torque) {
      t[i][0] = t_level[i][ilevel][0];
      t[i][1] = t_level[i][ilevel][1];
      t[i][2] = t_level[i][ilevel][2];
    }
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
  double ***t_level = fix_respa->t_level;
  double **t = atom->torque;
  int n = atom->nlocal;

  for (int ilevel = 1; ilevel < nlevels; ilevel++) {
    for (int i = 0; i < n; i++) {
      f[i][0] += f_level[i][ilevel][0];
      f[i][1] += f_level[i][ilevel][1];
      f[i][2] += f_level[i][ilevel][2];
      if (fix_respa->store_torque) {
        t[i][0] += t_level[i][ilevel][0];
        t[i][1] += t_level[i][ilevel][1];
        t[i][2] += t_level[i][ilevel][2];
      }
    }
  }
}

/*-----------------------------------------------------------------------
  set flags for when some hybrid forces should be computed
------------------------------------------------------------------------- */

void Respa::set_compute_flags(int ilevel)
{

  if (nhybrid_styles < 1) return;

  pair_compute = 0;
  for (int i=0; i<nhybrid_styles; ++i) {
    hybrid_compute[i] = (hybrid_level[i] == ilevel) ? 1 : 0;
    if (hybrid_compute[i]) pair_compute = 1;
  }
  tally_global = (ilevel == nlevels-1) ? 1 : 0;
}
