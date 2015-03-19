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
   Contributing authors: Paolo Raiteri (Curtin U) & Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_temp_stochastic.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "random_mars.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTempStochastic::FixTempStochastic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix temp/stochastic command");

  // Stochastic thermostat should be applied every step

  nevery = 1;
  scalar_flag = 1;
  global_freq = nevery;
  dynamic_group_allow = 1;
  extscalar = 1;

  tstr = NULL;
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[3][2]);
    tstyle = EQUAL;
  } else {
    t_start = force->numeric(FLERR,arg[3]);
    t_target = t_start;
    tstyle = CONSTANT;
  }

  t_stop = force->numeric(FLERR,arg[4]);
  t_period = force->numeric(FLERR,arg[5]);
  int seed = force->inumeric(FLERR,arg[6]);

  // error checks

  if (t_period <= 0.0) error->all(FLERR,"Fix temp/stochastic period must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix temp/stochastic random seed");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTempStochastic::~FixTempStochastic()
{
  delete [] tstr;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;

  delete random;
}

/* ---------------------------------------------------------------------- */

int FixTempStochastic::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempStochastic::init()
{

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix temp/stochastic does not exist");
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else error->all(FLERR,"Variable for fix temp/stochastic is invalid style");
  }

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/stochastic does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempStochastic::end_of_step()
{

  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  double delta = update->ntimestep - update->beginstep;

  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR,
                 "Fix temp/stochastic variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  t_current = temperature->compute_scalar();

  if (t_current == 0.0)
    error->all(FLERR,"Computed temperature for fix temp/stochastic cannot be 0.0");

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double tfactor = force->mvv2e / (3 * nlocal * force->boltz);
  double efactor = 0.5 * force->boltz * 3 * nlocal;

  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	  mass[type[i]];
  }

  t *= tfactor;

  if (t == 0.0)
    error->all(FLERR,"Computed temperature for fix temp/stochastic cannot be 0.0");

  //  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop-t_start);

  double ekin0 = efactor * t;
  double kbt = 0.5 * force->boltz * t_target;

  double lambda = resamplekin(kbt, ekin0, 3*nlocal, t_period/update->dt);

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	v[i][0] *= lambda;
	v[i][1] *= lambda;
	v[i][2] *= lambda;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	temperature->remove_bias(i,v[i]);
	v[i][0] *= lambda;
	v[i][1] *= lambda;
	v[i][2] *= lambda;
	temperature->restore_bias(i,v[i]);
      }
    }
  }

  if (group->count(igroup) == 0)
    error->all(FLERR,"Cannot zero momentum of 0 atoms");

  // compute velocity of center-of-mass of group

  double masstotal = group->mass(igroup);
  double vcm[3];
  group->vcm(igroup,masstotal,vcm);

  // adjust velocities by vcm to zero linear momentum

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] -= vcm[0];
      v[i][1] -= vcm[1];
      v[i][2] -= vcm[2];
    }

  double t_current = temperature->compute_scalar();
  energy += t_current * efactor;

}

/* ---------------------------------------------------------------------- */

int FixTempStochastic::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixTempStochastic::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

double FixTempStochastic::compute_scalar()
{
  return energy;
}

//////////////////////////////////////////////////////////////////////////

double FixTempStochastic::resamplekin(double kbt, double ekin0, int ndeg, double taut){
/*
  ndeg:  number of degrees of freedom of the atoms to be thermalized
  taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
*/
  double c1, c2, r1, r2, fscale2;
 
  if(taut>0.1){
    c1=exp(-1.0/taut);
  } else{
    c1=0.0;
  }
  c2 = (1.0-c1)*kbt/ekin0;
  r1 = random->gaussian();
  r2 = resamplekin_sumnoises(ndeg-1);
  fscale2 = c1 + c2*(r1*r1+r2) + 2.0*r1*sqrt(c1*c2);
  return sqrt(fscale2);

}

double FixTempStochastic::resamplekin_sumnoises(int nn){
/*
  returns the sum of n independent gaussian noises squared
   (i.e. equivalent to summing the square of the return values of nn calls to random->gaussian)
*/
  double rr;

  if(nn==0) {
    return 0.0;
  } else if(nn==1) {
    rr=random->gaussian();
    return rr*rr;
  } else if(nn%2==0) {
    return 2.0*gamdev(nn/2);
  } else {
    rr=random->gaussian();
    return 2.0*gamdev((nn-1)/2) + rr*rr;
  }
}

double FixTempStochastic::gamdev(const int ia)
{
  int j;
  double am,e,s,v1,v2,x,y;
  
  if (ia < 1) {}; // FATAL ERROR
  if (ia < 6) {
    x=1.0;
    for (j=1;j<=ia;j++) x *= random->uniform();
    x = -log(x);
  } else {
    restart:
    do {
      do {
      	do {
          v1=random->uniform();
          v2=2.0*random->uniform()-1.0;
      	} while (v1*v1+v2*v2 > 1.0);
      	y=v2/v1;
      	am=ia-1;
      	s=sqrt(2.0*am+1.0);
      	x=s*y+am;
      } while (x <= 0.0);
      if (am*log(x/am)-s*y < -700 || v1<0.00001) {
        goto restart;
      }
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    } while (random->uniform() > e);
  }
  return x;
}

