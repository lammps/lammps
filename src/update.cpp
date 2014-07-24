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

#include "string.h"
#include "stdlib.h"
#include "update.h"
#include "integrate.h"
#include "min.h"
#include "style_integrate.h"
#include "style_minimize.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "region.h"
#include "compute.h"
#include "output.h"
#include "memory.h"
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
  setupflag = 0;
  multireplica = 0;

  restrict_output = 0;

  eflag_global = vflag_global = -1;

  unit_style = NULL;
  set_units("lj");

  integrate_style = NULL;
  integrate = NULL;
  minimize_style = NULL;
  minimize = NULL;

  str = (char *) "verlet";
  create_integrate(1,&str,lmp->suffix);

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
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  // if USER-CUDA mode is enabled:
  // integrate/minimize style must be CUDA variant

  if (whichflag == 1 && lmp->cuda)
    if (strstr(integrate_style,"cuda") == NULL)
      error->all(FLERR,"USER-CUDA mode requires CUDA variant of run style");
  if (whichflag == 2 && lmp->cuda)
    if (strstr(minimize_style,"cuda") == NULL)
      error->all(FLERR,"USER-CUDA mode requires CUDA variant of min style");

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
    force->hplanck = 0.18292026;  // using LJ parameters for argon
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
    force->mv2d = 1.0 / 0.602214179;
    force->nktv2p = 68568.415;
    force->qqr2e = 332.06371;
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
    force->mv2d = 1.0 / 0.602214179;
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
    force->femtosecond = 0.0241888428;
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

void Update::create_integrate(int narg, char **arg, char *suffix)
{
  if (narg < 1) error->all(FLERR,"Illegal run_style command");

  delete [] integrate_style;
  delete integrate;

  int sflag;
  new_integrate(arg[0],narg-1,&arg[1],suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",arg[0],suffix);
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
                           char *suffix, int &sflag)
{
  int success = 0;

  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    success = 1;

    if (0) return;

#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
    else if (strcmp(estyle,#key) == 0) integrate = new Class(lmp,narg,arg);
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS

    else success = 0;
  }

  if (!success) {
    sflag = 0;

    if (0) return;

#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
    else if (strcmp(style,#key) == 0) integrate = new Class(lmp,narg,arg);
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS

    else error->all(FLERR,"Illegal integrate style");
  }
}

/* ---------------------------------------------------------------------- */

void Update::create_minimize(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal min_style command");

  delete [] minimize_style;
  delete minimize;

  if (0) return;      // dummy line to enable else-if macro expansion

#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) minimize = new Class(lmp);
#include "style_minimize.h"
#undef MINIMIZE_CLASS

  else error->all(FLERR,"Illegal min_style command");

  int n = strlen(arg[0]) + 1;
  minimize_style = new char[n];
  strcpy(minimize_style,arg[0]);
}

/* ----------------------------------------------------------------------
   reset timestep as called from input script
------------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal reset_timestep command");
  bigint newstep = ATOBIGINT(arg[0]);
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
  if (ntimestep > MAXBIGINT) error->all(FLERR,"Too big a timestep");

  // set atimestep to new timestep
  // so future update_time() calls will be correct

  atimestep = ntimestep;

  // trigger reset of timestep for output and for fixes that require it
  // do not allow any timestep-dependent fixes to be defined

  output->reset_timestep(ntimestep);

  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->time_depend)
      error->all(FLERR,
                 "Cannot reset timestep with a time-dependent fix defined");
    modify->fix[i]->reset_timestep(ntimestep);
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

  // set last_build of all neigh lists to -1 to force rebuild

  for (int i = 0; i < neighbor->nlist; i++)
    neighbor->lists[i]->last_build = -1;

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
