/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "update.h"

#include "style_integrate.h"    // IWYU pragma: keep
#include "style_minimize.h"     // IWYU pragma: keep

#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "integrate.h"
#include "min.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"

#include <cstring>

using namespace LAMMPS_NS;

// template for factory functions:
// there will be one instance for each style keyword in the respective style_xxx.h files

template <typename T> static Integrate *integrate_creator(LAMMPS *lmp, int narg, char **arg)
{
  return new T(lmp, narg, arg);
}

template <typename T> static Min *minimize_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ---------------------------------------------------------------------- */

Update::Update(LAMMPS *lmp) :
    Pointers(lmp), unit_style(nullptr), integrate(nullptr), integrate_style(nullptr),
    minimize(nullptr), minimize_style(nullptr), integrate_map(nullptr), minimize_map(nullptr)
{
  char *str;

  ntimestep = 0;
  atime = 0.0;
  atimestep = 0;
  first_update = 0;

  whichflag = 0;
  firststep = laststep = 0;
  beginstep = endstep = 0;
  restrict_output = 0;
  setupflag = 0;
  multireplica = 0;

  eflag_global = vflag_global = -1;
  eflag_atom = vflag_atom = 0;

  dt_default = 1;
  dt = 0.0;
  set_units("lj");

  integrate_map = new IntegrateCreatorMap();

#define INTEGRATE_CLASS
#define IntegrateStyle(key, Class) (*integrate_map)[#key] = &integrate_creator<Class>;
#include "style_integrate.h"    // IWYU pragma: keep
#undef IntegrateStyle
#undef INTEGRATE_CLASS

  minimize_map = new MinimizeCreatorMap();

#define MINIMIZE_CLASS
#define MinimizeStyle(key, Class) (*minimize_map)[#key] = &minimize_creator<Class>;
#include "style_minimize.h"    // IWYU pragma: keep
#undef MinimizeStyle
#undef MINIMIZE_CLASS

  str = (char *) "verlet";
  create_integrate(1, &str, 1);

  str = (char *) "cg";
  create_minimize(1, &str, 1);
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  delete[] unit_style;

  delete[] integrate_style;
  delete integrate;

  delete[] minimize_style;
  delete minimize;

  delete integrate_map;
  delete minimize_map;
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  // init the appropriate integrate and/or minimize class
  // if neither (e.g. from write_restart) then just return

  if (whichflag == 0) return;
  if (whichflag == 1)
    integrate->init();
  else if (whichflag == 2)
    minimize->init();

  // only set first_update if a run or minimize is being performed

  first_update = 1;
}

/* ---------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
  // physical constants from:
  // https://physics.nist.gov/cuu/Constants/Table/allascii.txt
  // using thermochemical calorie = 4.184 J

  double dt_old = dt;

  if (strcmp(style, "lj") == 0) {
    force->boltz = 1.0;
    force->hplanck = 1.0;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;

    dt = 0.005;
    neighbor->skin = 0.3;

  } else if (strcmp(style, "real") == 0) {
    force->boltz = 0.0019872067;
    force->hplanck = 95.306976368;
    force->mvv2e = 48.88821291 * 48.88821291;
    force->ftm2v = 1.0 / 48.88821291 / 48.88821291;
    force->mv2d = 1.0 / 0.602214129;
    force->nktv2p = 68568.415;
    force->qqr2e = 332.06371;    // see also force->qqr2d_lammps_real
    force->qe2f = 23.060549;
    force->vxmu2f = 1.4393264316e4;
    force->xxt2kmu = 0.1;
    force->e_mass = 1.0 / 1836.1527556560675;
    force->hhmrr2e = 0.0957018663603261;
    force->mvh2r = 1.5339009481951;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;

    dt = 1.0;
    neighbor->skin = 2.0;

  } else if (strcmp(style, "metal") == 0) {
    force->boltz = 8.617343e-5;
    force->hplanck = 4.135667403e-3;
    force->mvv2e = 1.0364269e-4;
    force->ftm2v = 1.0 / 1.0364269e-4;
    force->mv2d = 1.0 / 0.602214129;
    force->nktv2p = 1.6021765e6;
    force->qqr2e = 14.399645;
    force->qe2f = 1.0;
    force->vxmu2f = 0.6241509647;
    force->xxt2kmu = 1.0e-4;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0e-3;
    force->qelectron = 1.0;

    dt = 0.001;
    neighbor->skin = 2.0;

  } else if (strcmp(style, "si") == 0) {
    force->boltz = 1.3806504e-23;
    force->hplanck = 6.62606896e-34;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.9876e9;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-10;
    force->femtosecond = 1.0e-15;
    force->qelectron = 1.6021765e-19;

    dt = 1.0e-8;
    neighbor->skin = 0.001;

  } else if (strcmp(style, "cgs") == 0) {
    force->boltz = 1.3806504e-16;
    force->hplanck = 6.62606896e-27;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-8;
    force->femtosecond = 1.0e-15;
    force->qelectron = 4.8032044e-10;

    dt = 1.0e-8;
    neighbor->skin = 0.1;

  } else if (strcmp(style, "electron") == 0) {
    force->boltz = 3.16681534e-6;
    force->hplanck = 0.1519829846;
    force->mvv2e = 1.06657236;
    force->ftm2v = 0.937582899;
    force->mv2d = 1.0;
    force->nktv2p = 2.94210108e13;
    force->qqr2e = 1.0;
    force->qe2f = 1.94469051e-10;
    force->vxmu2f = 3.39893149e1;
    force->xxt2kmu = 3.13796367e-2;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.88972612;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;

    dt = 0.001;
    neighbor->skin = 2.0;

  } else if (strcmp(style, "micro") == 0) {
    force->boltz = 1.3806504e-8;
    force->hplanck = 6.62606896e-13;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.987556e6;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-4;
    force->femtosecond = 1.0e-9;
    force->qelectron = 1.6021765e-7;

    dt = 2.0;
    neighbor->skin = 0.1;

  } else if (strcmp(style, "nano") == 0) {
    force->boltz = 0.013806504;
    force->hplanck = 6.62606896e-4;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 230.7078669;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-1;
    force->femtosecond = 1.0e-6;
    force->qelectron = 1.0;

    dt = 0.00045;
    neighbor->skin = 0.1;

  } else
    error->all(FLERR, "Illegal units command");

  delete[] unit_style;
  unit_style = utils::strdup(style);

  // check if timestep was changed from default value
  if (!dt_default && (comm->me == 0)) {
    error->warning(FLERR, "Changing timestep from {:.6} to {:.6} due to changing units to {}",
                   dt_old, dt, unit_style);
  }
  dt_default = 1;
}

/* ---------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg, int trysuffix)
{
  if (narg < 1) error->all(FLERR, "Illegal run_style command");

  delete[] integrate_style;
  delete integrate;
  integrate_style = nullptr;
  integrate = nullptr;

  int sflag;

  if (narg - 1 > 0) {
    new_integrate(arg[0], narg - 1, &arg[1], trysuffix, sflag);
  } else {
    new_integrate(arg[0], 0, nullptr, trysuffix, sflag);
  }

  std::string estyle = arg[0];
  if (sflag) {
    estyle += "/";
    if (sflag == 1)
      estyle += lmp->suffix;
    else if (sflag == 2)
      estyle += lmp->suffix2;
    else if ((sflag == 3) && lmp->non_pair_suffix())
      estyle += lmp->non_pair_suffix();
  }
  integrate_style = utils::strdup(estyle);
}

/* ----------------------------------------------------------------------
   create the Integrate style, first with suffix appended
------------------------------------------------------------------------- */

void Update::new_integrate(char *style, int narg, char **arg, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2 * lmp->pair_only_flag;
      std::string estyle = style + std::string("/") + lmp->non_pair_suffix();
      if (integrate_map->find(estyle) != integrate_map->end()) {
        IntegrateCreator &integrate_creator = (*integrate_map)[estyle];
        integrate = integrate_creator(lmp, narg, arg);
        return;
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + std::string("/") + lmp->suffix2;
      if (integrate_map->find(estyle) != integrate_map->end()) {
        IntegrateCreator &integrate_creator = (*integrate_map)[estyle];
        integrate = integrate_creator(lmp, narg, arg);
        return;
      }
    }
  }

  sflag = 0;
  if (integrate_map->find(style) != integrate_map->end()) {
    IntegrateCreator &integrate_creator = (*integrate_map)[style];
    integrate = integrate_creator(lmp, narg, arg);
    return;
  }

  error->all(FLERR, "Illegal integrate style");
}

/* ---------------------------------------------------------------------- */

void Update::create_minimize(int narg, char **arg, int trysuffix)
{
  if (narg < 1) error->all(FLERR, "Illegal minimize_style command");

  delete[] minimize_style;
  delete minimize;
  minimize_style = nullptr;
  minimize = nullptr;

  // temporarily assign the style name without suffix (for error messages during creation)
  minimize_style = arg[0];

  int sflag;
  new_minimize(arg[0], narg - 1, &arg[1], trysuffix, sflag);

  std::string estyle = arg[0];
  if (sflag) {
    estyle += "/";
    if (sflag == 1)
      estyle += lmp->suffix;
    else if (sflag == 2)
      estyle += lmp->suffix2;
    else if ((sflag == 3) && lmp->non_pair_suffix())
      estyle += lmp->non_pair_suffix();
  }
  minimize_style = utils::strdup(estyle);
}

/* ----------------------------------------------------------------------
   create the Minimize style, first with suffix appended
------------------------------------------------------------------------- */

void Update::new_minimize(char *style, int /* narg */, char ** /* arg */, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2 * lmp->pair_only_flag;
      std::string estyle = style + std::string("/") + lmp->non_pair_suffix();
      if (minimize_map->find(estyle) != minimize_map->end()) {
        MinimizeCreator &minimize_creator = (*minimize_map)[estyle];
        minimize = minimize_creator(lmp);
        return;
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + std::string("/") + lmp->suffix2;
      if (minimize_map->find(estyle) != minimize_map->end()) {
        MinimizeCreator &minimize_creator = (*minimize_map)[estyle];
        minimize = minimize_creator(lmp);
        return;
      }
    }
  }

  sflag = 0;
  if (minimize_map->find(style) != minimize_map->end()) {
    MinimizeCreator &minimize_creator = (*minimize_map)[style];
    minimize = minimize_creator(lmp);
    return;
  }

  error->all(FLERR, "Illegal minimize style");
}

/* ----------------------------------------------------------------------
   reset timestep called from input script
------------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "reset_timestep", error);

  reset_timestep(utils::bnumeric(FLERR, arg[0], false, lmp), true);

  if (narg > 1) {
    if (strcmp(arg[1], "time") == 0) {
      if (narg < 3) utils::missing_cmd_args(FLERR, "reset_timestep time", error);
      atimestep = ntimestep;
      atime = utils::numeric(FLERR, arg[2], false, lmp);
    } else
      error->all(FLERR, "Unknown reset_timestep option {}", arg[1]);
  }
}

/* ----------------------------------------------------------------------
   reset timestep
   called from input script (indirectly) or rerun command
------------------------------------------------------------------------- */

void Update::reset_timestep(bigint newstep, bool do_check)
{
  if (newstep < 0) error->all(FLERR, "Timestep must be >= 0");

  bigint oldstep = ntimestep;
  ntimestep = newstep;

  // if newstep >= oldstep, update simulation time accordingly
  // if newstep < oldstep, zero simulation time

  if (newstep >= oldstep) update_time();

  if (newstep < oldstep) {
    atime = 0.0;
    atimestep = newstep;
  }

  // changes to output that depend on timestep
  // no active dumps allowed

  output->reset_timestep(ntimestep);

  // rerun will not be meaningful with this check active.
  if (do_check) {
    // do not allow timestep-dependent fixes to be defined

    for (const auto &ifix : modify->get_fix_list())
      if (ifix->time_depend)
        error->all(FLERR, "Cannot reset timestep with time-dependent fix {} defined", ifix->style);
  }

  // reset eflag/vflag global so no commands will think eng/virial are current

  eflag_global = vflag_global = -1;

  // reset invoked flags of computes, so no commands will think they are current between runs
  // clear timestep list of computes that store future invocation times

  for (const auto &icompute : modify->get_compute_list()) {
    icompute->invoked_scalar = -1;
    icompute->invoked_vector = -1;
    icompute->invoked_array = -1;
    icompute->invoked_peratom = -1;
    icompute->invoked_local = -1;
    if (icompute->timeflag) icompute->clearstep();
  }

  // neighbor Bin/Stencil/Pair classes store timestamps that need to be cleared

  neighbor->reset_timestep(ntimestep);
}

/* ----------------------------------------------------------------------
   update elapsed simulation time
   called at end of runs or when timestep size changes
------------------------------------------------------------------------- */

void Update::update_time()
{
  atime += (ntimestep - atimestep) * dt;
  atimestep = ntimestep;
}

/* ----------------------------------------------------------------------
   memory usage of update and integrate/minimize
------------------------------------------------------------------------- */

double Update::memory_usage()
{
  double bytes = 0;
  if (whichflag == 1)
    bytes += integrate->memory_usage();
  else if (whichflag == 2)
    bytes += minimize->memory_usage();
  return bytes;
}
