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
#include "style_integrate.h"
#include "style_minimize.h"
#include "neighbor.h"
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
  int n;
  char *str;

  ntimestep = 0;
  first_update = 0;

  whichflag = 0;
  firststep = laststep = 0;
  beginstep = endstep = 0;
  setupflag = 0;

  restrict_output = 0;

  eflag_global = vflag_global = -1;

  unit_style = NULL;
  set_units("lj");

  str = (char *) "verlet";
  n = strlen(str) + 1;
  integrate_style = new char[n];
  strcpy(integrate_style,str);
  integrate = new Verlet(lmp,0,NULL);

  str = (char *) "cg";
  n = strlen(str) + 1;
  minimize_style = new char[n];
  strcpy(minimize_style,str);
  minimize = new MinCG(lmp);
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
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    dt = 0.005;
    neighbor->skin = 0.3;
    
  } else if (strcmp(style,"real") == 0) {
    force->boltz = 0.0019872067;
    force->mvv2e = 48.88821291 * 48.88821291;
    force->ftm2v = 1.0 / 48.88821291 / 48.88821291;
    force->mv2d = 1.0 / 0.602214179;
    force->nktv2p = 68568.415;
    force->qqr2e = 332.06371;
    force->qe2f = 23.060549; 
    force->vxmu2f = 1.4393264316e4;
    force->xxt2kmu = 0.1;
    dt = 1.0;
    neighbor->skin = 2.0;

  } else if (strcmp(style,"metal") == 0) {
    force->boltz = 8.617343e-5;
    force->mvv2e = 1.0364269e-4;
    force->ftm2v = 1.0 / 1.0364269e-4;
    force->mv2d = 1.0 / 0.602214179;
    force->nktv2p = 1.6021765e6;
    force->qqr2e = 14.399645;
    force->qe2f = 1.0;
    force->vxmu2f = 0.6241509647;
    force->xxt2kmu = 1.0e-4;
    dt = 0.001;
    neighbor->skin = 2.0;

  } else if (strcmp(style,"si") == 0) {
    force->boltz = 1.3806504e-23;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.9876e9;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    dt = 1.0e-8;
    neighbor->skin = 0.001;

  } else if (strcmp(style,"cgs") == 0) {
    force->boltz = 1.3806504e-16;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    dt = 1.0e-8;
    neighbor->skin = 0.1;

  } else if (strcmp(style,"electron") == 0) {
    force->boltz = 3.16679e-6;
    force->mvv2e = 1.0;
    force->ftm2v = 1 / 1.0327499;	
    force->mv2d = 1.0;
    force->nktv2p = 2.9421e13;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 0.6241509647;
    force->xxt2kmu = 1.0e-4;
    dt = 0.001;
    neighbor->skin = 2.0;

  } else error->all("Illegal units command");

  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal run_style command");

  delete [] integrate_style;
  delete integrate;

  if (0) return;      // dummy line to enable else-if macro expansion

#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) integrate = new Class(lmp,narg-1,&arg[1]);
#include "style_integrate.h"
#undef INTEGRATE_CLASS

  else error->all("Illegal run_style command");

  int n = strlen(arg[0]) + 1;
  integrate_style = new char[n];
  strcpy(integrate_style,arg[0]);
}

/* ---------------------------------------------------------------------- */

void Update::create_minimize(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal min_style command");

  delete [] minimize_style;
  delete minimize;

  if (0) return;      // dummy line to enable else-if macro expansion

#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) minimize = new Class(lmp);
#include "style_minimize.h"
#undef MINIMIZE_CLASS

  else error->all("Illegal min_style command");

  int n = strlen(arg[0]) + 1;
  minimize_style = new char[n];
  strcpy(minimize_style,arg[0]);
}

/* ----------------------------------------------------------------------
   reset timestep from input script
   do not allow dump files or a restart to be defined
   do not allow any timestep-dependent fixes to be defined
   do not allow any dynamic regions to be defined
   reset eflag/vflag global so nothing will think eng/virial are current
   reset invoked flags of computes,
     so nothing will think they are current between runs
   clear timestep list of computes that store future invocation times
------------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal reset_timestep command");

  for (int i = 0; i < output->ndump; i++)
    if (output->last_dump[i] >= 0)
      error->all("Cannot reset timestep with dump file already written to");
  if (output->restart && output->last_restart >= 0)
    error->all("Cannot reset timestep with restart file already written");

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->time_depend)
      error->all("Cannot reset timestep with a time-dependent fix defined");

  for (int i = 0; i < domain->nregion; i++)
    if (domain->regions[i]->dynamic_check())
      error->all("Cannot reset timestep with a dynamic region defined");

  eflag_global = vflag_global = -1;

  for (int i = 0; i < modify->ncompute; i++) {
    modify->compute[i]->invoked_scalar = -1;
    modify->compute[i]->invoked_vector = -1;
    modify->compute[i]->invoked_array = -1;
    modify->compute[i]->invoked_peratom = -1;
    modify->compute[i]->invoked_local = -1;
  }

  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();

  ntimestep = atoi(arg[0]);
  if (ntimestep < 0) error->all("Timestep must be >= 0");
}

/* ----------------------------------------------------------------------
   memory usage of update and integrate/minimize
------------------------------------------------------------------------- */

double Update::memory_usage()
{
  double bytes = 0.0;
  if (whichflag == 1) bytes += integrate->memory_usage();
  else if (whichflag == 2) bytes += minimize->memory_usage();
  return bytes;
}
