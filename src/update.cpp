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

#include "update.h"
#include <cstring>
#include "integrate.h"
#include "min.h"
#include "style_integrate.h"
#include "style_minimize.h"
#include "neighbor.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Update::Update(LAMMPS *lmp) : Pointers(lmp)
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
  post_integrate = 0;
  multireplica = 0;

  eflag_global = vflag_global = -1;

  unit_style = NULL;
  set_units("lj");

  integrate_style = NULL;
  integrate = NULL;
  minimize_style = NULL;
  minimize = NULL;

  integrate_map = new IntegrateCreatorMap();

#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
  (*integrate_map)[#key] = &integrate_creator<Class>;
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS

  minimize_map = new MinimizeCreatorMap();

#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) \
  (*minimize_map)[#key] = &minimize_creator<Class>;
#include "style_minimize.h"
#undef MinimizeStyle
#undef MINIMIZE_CLASS

  str = (char *) "verlet";
  create_integrate(1,&str,1);

  str = (char *) "cg";
  create_minimize(1,&str);
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  delete [] unit_style;

  delete [] integrate_style;
  delete integrate;

  delete [] minimize_style;
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
  if (whichflag == 1) integrate->init();
  else if (whichflag == 2) minimize->init();

  // only set first_update if a run or minimize is being performed

  first_update = 1;
}

/* ---------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
  // physical constants from:
  // http://physics.nist.gov/cuu/Constants/Table/allascii.txt
  // using thermochemical calorie = 4.184 J

  if (strcmp(style,"lj") == 0) {
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

  } else if (strcmp(style,"real") == 0) {
    force->boltz = 0.0019872067;
    force->hplanck = 95.306976368;
    force->mvv2e = 48.88821291 * 48.88821291;
    force->ftm2v = 1.0 / 48.88821291 / 48.88821291;
    force->mv2d = 1.0 / 0.602214129;
    force->nktv2p = 68568.415;
    force->qqr2e = 332.06371;     // see also force->qqr2d_lammps_real
    force->qe2f = 23.060549;
    force->vxmu2f = 1.4393264316e4;
    force->xxt2kmu = 0.1;
    force->e_mass = 1.0/1836.1527556560675;
    force->hhmrr2e = 0.0957018663603261;
    force->mvh2r = 1.5339009481951;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;

    dt = 1.0;
    neighbor->skin = 2.0;

  } else if (strcmp(style,"metal") == 0) {
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

  } else if (strcmp(style,"si") == 0) {
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

  } else if (strcmp(style,"cgs") == 0) {
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

  } else if (strcmp(style,"electron") == 0) {
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

  } else if (strcmp(style,"micro") == 0) {
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

  } else if (strcmp(style,"nano") == 0) {
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

  } else error->all(FLERR,"Illegal units command");

  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg, int trysuffix)
{
  if (narg < 1) error->all(FLERR,"Illegal run_style command");

  delete [] integrate_style;
  delete integrate;

  int sflag;
  new_integrate(arg[0],narg-1,&arg[1],trysuffix,sflag);

  if (sflag) {
    char estyle[256];
    if (sflag == 1) snprintf(estyle,256,"%s/%s",arg[0],lmp->suffix);
    else snprintf(estyle,256,"%s/%s",arg[0],lmp->suffix2);
    int n = strlen(estyle) + 1;
    integrate_style = new char[n];
    strcpy(integrate_style,estyle);
  } else {
    int n = strlen(arg[0]) + 1;
    integrate_style = new char[n];
    strcpy(integrate_style,arg[0]);
  }
}

/* ----------------------------------------------------------------------
   create the Integrate style, first with suffix appended
------------------------------------------------------------------------- */

void Update::new_integrate(char *style, int narg, char **arg,
                           int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      char estyle[256];
      snprintf(estyle,256,"%s/%s",style,lmp->suffix);
      if (integrate_map->find(estyle) != integrate_map->end()) {
        IntegrateCreator integrate_creator = (*integrate_map)[estyle];
        integrate = integrate_creator(lmp, narg, arg);
        return;
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      char estyle[256];
      snprintf(estyle,256,"%s/%s",style,lmp->suffix2);
      if (integrate_map->find(estyle) != integrate_map->end()) {
        IntegrateCreator integrate_creator = (*integrate_map)[estyle];
        integrate = integrate_creator(lmp, narg, arg);
        return;
      }
    }
  }

  sflag = 0;
  if (integrate_map->find(style) != integrate_map->end()) {
    IntegrateCreator integrate_creator = (*integrate_map)[style];
    integrate = integrate_creator(lmp, narg, arg);
    return;
  }

  error->all(FLERR,"Illegal integrate style");
}

/* ----------------------------------------------------------------------
   one instance per integrate style in style_integrate.h
------------------------------------------------------------------------- */

template <typename T>
Integrate *Update::integrate_creator(LAMMPS *lmp, int narg, char ** arg)
{
  return new T(lmp, narg, arg);
}

/* ---------------------------------------------------------------------- */

void Update::create_minimize(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal min_style command");

  delete [] minimize_style;
  delete minimize;

  if (minimize_map->find(arg[0]) != minimize_map->end()) {
    MinimizeCreator minimize_creator = (*minimize_map)[arg[0]];
    minimize = minimize_creator(lmp);
  }
  else error->all(FLERR,"Illegal min_style command");

  int n = strlen(arg[0]) + 1;
  minimize_style = new char[n];
  strcpy(minimize_style,arg[0]);
}

/* ----------------------------------------------------------------------
   one instance per minimize style in style_minimize.h
------------------------------------------------------------------------- */

template <typename T>
Min *Update::minimize_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   reset timestep as called from input script
------------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal reset_timestep command");
  bigint newstep = force->bnumeric(FLERR,arg[0]);
  reset_timestep(newstep);
}

/* ----------------------------------------------------------------------
   reset timestep
   called from rerun command and input script (indirectly)
------------------------------------------------------------------------- */

void Update::reset_timestep(bigint newstep)
{
  ntimestep = newstep;
  if (ntimestep < 0) error->all(FLERR,"Timestep must be >= 0");

  // set atimestep to new timestep
  // so future update_time() calls will be correct

  atimestep = ntimestep;

  // trigger reset of timestep for output
  // do not allow any timestep-dependent fixes to be already defined

  output->reset_timestep(ntimestep);

  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->time_depend)
      error->all(FLERR,
                 "Cannot reset timestep with a time-dependent fix defined");
  }

  // reset eflag/vflag global so no commands will think eng/virial are current

  eflag_global = vflag_global = -1;

  // reset invoked flags of computes,
  // so no commands will think they are current between runs

  for (int i = 0; i < modify->ncompute; i++) {
    modify->compute[i]->invoked_scalar = -1;
    modify->compute[i]->invoked_vector = -1;
    modify->compute[i]->invoked_array = -1;
    modify->compute[i]->invoked_peratom = -1;
    modify->compute[i]->invoked_local = -1;
  }

  // clear timestep list of computes that store future invocation times

  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();

  // Neighbor Bin/Stencil/Pair classes store timestamps that need to be cleared

  neighbor->reset_timestep(ntimestep);

  // NOTE: 7Jun12, adding rerun command, don't think this is required

  //for (int i = 0; i < domain->nregion; i++)
  //  if (domain->regions[i]->dynamic_check())
  //    error->all(FLERR,"Cannot reset timestep with a dynamic region defined");
}

/* ----------------------------------------------------------------------
   update elapsed simulation time
   called at end of runs or when timestep size changes
------------------------------------------------------------------------- */

void Update::update_time()
{
  atime += (ntimestep-atimestep) * dt;
  atimestep = ntimestep;
}

/* ----------------------------------------------------------------------
   memory usage of update and integrate/minimize
------------------------------------------------------------------------- */

bigint Update::memory_usage()
{
  bigint bytes = 0;
  if (whichflag == 1) bytes += integrate->memory_usage();
  else if (whichflag == 2) bytes += minimize->memory_usage();
  return bytes;
}
